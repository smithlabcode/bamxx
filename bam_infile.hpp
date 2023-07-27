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

#ifndef BAM_INFILE_HPP
#define BAM_INFILE_HPP

#include <string>

#include <htslib/sam.h>

#include "bam_record.hpp"

class bam_infile {
public:
  bam_infile(): fp(nullptr), error_code(0) {}

  ~bam_infile() {
    if (fp != nullptr) hts_close(fp);
    if (hdr != nullptr) sam_hdr_destroy(hdr);
  }

  bam_infile(const std::string &filename) {
    fp = hts_open(filename.c_str(), "r");
    error_code = (fp) ? 0 : -1;
    if (!error_code) {
      hdr = sam_hdr_read(fp);
      if (!hdr) {
        if (fp != nullptr) hts_close(fp);
        error_code = -2;
      }
    }
    else if (fp != nullptr)
      hts_close(fp);
  }

  const sam_hdr_t *const get_hdr() const { return hdr; };
  //MN: temporary solution
  sam_hdr_t * get_hdr_nonconst() { return hdr; };

  samFile *const get_fp() const { return fp; };
  samFile *get_ptr() { return fp; };

  bam_infile &get_bam_rec(bam_rec &br) {
    if (error_code) return *this; // already in error; do nothing

    bam1_t *aln = bam_init1();
    error_code = sam_read1(fp, hdr, aln);
    if (error_code < 0) {
      bam_destroy1(aln);
      if (error_code < -1)
        throw std::runtime_error("failure reading BAM/SAM file");
    }
    else {
      error_code = 0; // ADS: `sam_read1` returns positive??
      bam_rec::steal(aln, br);
    }
    return *this;
  }
  bool is_mapped_reads_file() const {
    const htsFormat *fmt = hts_get_format(fp);
    return fmt->category == sequence_data &&
      (fmt->format == bam || fmt->format == sam);
  }

  operator bool() const { return (error_code == 0); }
private:
  samFile *fp{};
  sam_hdr_t *hdr{};
  int error_code{};
};

inline bam_infile &
get_bam(bam_infile &in, bam_rec &br) {
  return in.get_bam_rec(br);
}

#endif
