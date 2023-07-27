/* MIT License
 *
 * Copyright (c) 2023 Andrew Smith and Masaru Nakajima
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef BAM_OUTFILE_HPP
#define BAM_OUTFILE_HPP

#include <string>
#include <stdexcept>

#include <htslib/sam.h>

#include "bam_record.hpp"

class bam_outfile {
public:
  bam_outfile(): fp(nullptr), hdr(nullptr), error_code(0) {}

  bam_outfile(const std::string &filename, const bam_header &init_hdr,
              const bool bam_fmt = false) {
    if (init_hdr) {
      fp = hts_open(filename.c_str(), (bam_fmt ? "bw" : "w"));
      error_code = (fp) ? 0 : -1;
      if (!error_code) {
        hdr = sam_hdr_dup(init_hdr.h);
        if (!hdr)
          error_code = -2;
        else {
          error_code = sam_hdr_write(fp, hdr);
          if (error_code > 0) error_code = 0;
        }
      }
    }
    else
      error_code = -2;
  }

  ~bam_outfile() {
    if (fp != nullptr) {
      hts_close(fp);
      fp = nullptr;
    }
    if (hdr != nullptr) {
      bam_hdr_destroy(hdr);
      hdr = nullptr;
    }
  }

  samFile *const get_fp() const { return fp; };
  samFile *get_ptr() { return fp; };

  bam_outfile &put_bam_rec(const bam_rec &br) {
    if (!error_code) { // already in error; do nothing
      error_code = sam_write1(fp, hdr, &br.get_rec());
      if (error_code > 0) error_code = 0;
      if (error_code < 0)
        throw std::runtime_error("failed writing to BAM/SAM file");
    }
    return *this;
  }

  operator bool() const { return (error_code == 0); }

private:
  htsFile *fp{};
  sam_hdr_t *hdr{};
  int error_code{};
};

bam_outfile &
operator<<(bam_outfile &out, const bam_rec &br) {
  return out.put_bam_rec(br);
}

#endif
