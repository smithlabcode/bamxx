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

  samFile *const get_fp() const { return fp; };

  bam_infile &get_bam_rec(bam_rec &br) {
    if (error_code) return *this; // already in error; do nothing

    bam1_t *aln = bam_init1();
    error_code = sam_read1(fp, hdr, aln);
    if (error_code < 0)
      bam_destroy1(aln);
    else {
      error_code = 0; // ADS: `sam_read1` returns positive??
      bam_rec::steal(aln, br);
    }
    return *this;
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
