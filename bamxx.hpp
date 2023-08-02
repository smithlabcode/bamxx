/* MIT License
 *
 * Copyright (c) 2023 Andrew D Smith and Masaru Nakajima
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* ADS: need to control all the macros from HTSlib pollution. For
   functions maybe:

   $ gcc -dM -E sam.h | grep "define [a-z]" | awk '{print $2}' |\
       grep "[(]" | awk -v FS="(" '{print "#undef",$1}'

   This gives about 65 symbols that need to be deleted. For the others
   I don't know what to do because some of them have "#define _" which
   means they should be system symbols.
*/

#ifndef BAM_RECORD_HPP
#define BAM_RECORD_HPP

#include <htslib/sam.h>          // from HTSlib
#include <htslib/thread_pool.h>  // from HTSlib

#include <stdexcept>
#include <string>
#include <utility>

namespace bamxx {

  struct bam_rec {
    bam_rec() = default;

    bam_rec(const bam_rec &other): b{bam_copy1(bam_init1(), other.b)} {}

    auto operator=(bam_rec rhs) -> bam_rec & {
      std::swap(b, rhs.b);
      return *this;
    }

    ~bam_rec() {
      if (b != nullptr) bam_destroy1(b);
    }

    bam1_t *b{};
  };

  struct bam_infile {
    explicit bam_infile(const std::string &fn): f{hts_open(fn.c_str(), "r")} {}

    ~bam_infile() {
      if (f != nullptr) hts_close(f);
    }

    operator bool() const { return f != nullptr; }

    template<typename T> auto read(T &h, bam_rec &b) -> bool {
      const int x = sam_read1(f, h.h, b.b);  // -1 on EOF; args non-const
      if (x < -1) throw std::runtime_error("failed reading bam record");
      return x >= 0;
    }

    auto is_mapped_reads_file() const -> bool {
      const htsFormat *fmt = hts_get_format(f);
      return fmt->category == sequence_data &&
        (fmt->format == bam || fmt->format == sam);
    }

    samFile *f{};
  };

  struct bam_header {
    bam_header() = default;

    bam_header(const bam_header &rhs): h{bam_hdr_dup(rhs.h)} {}

    explicit bam_header(bam_infile &in): h{sam_hdr_read(in.f)} {}

    ~bam_header() {
      if (h != nullptr) bam_hdr_destroy(h);
    }

    operator bool() const { return h != nullptr; }

    auto add_pg_line(const std::string &cmd, const std::string &id,
                     const std::string &vn) -> bool {
      return sam_hdr_add_line(h, "PG", "ID", id.c_str(), "VN", vn.c_str(), "CL",
                              cmd.c_str(), nullptr) == 0;
    }

    auto tostring() const -> std::string { return sam_hdr_str(h); }

    sam_hdr_t *h{};
  };

  struct bam_outfile {
    explicit bam_outfile(const std::string &fn, const bool fmt = false)
      : f{hts_open(fn.c_str(), fmt ? "bw" : "w")} {}

    ~bam_outfile() {
      if (f != nullptr) hts_close(f);
    }

    operator bool() const { return f != nullptr; }

    auto write(const bam_header &h, const bam_rec &b) -> bool {
      return sam_write1(f, h.h, b.b) >= 0;
    }

    auto write(const bam_header &h) -> bool { return sam_hdr_write(f, h.h) == 0; }

    htsFile *f{};
  };

  struct bam_tpool {
    explicit bam_tpool(const int n_threads)
      : tpool{hts_tpool_init(n_threads), 0} {}

    ~bam_tpool() { hts_tpool_destroy(tpool.pool); }

    template<class T> auto set_io(const T &bam_file) -> void {
      const int ret = hts_set_thread_pool(bam_file.f, &tpool);
      if (ret < 0) throw std::runtime_error("failed to set thread pool");
    }

    htsThreadPool tpool{};
  };
};

#endif
