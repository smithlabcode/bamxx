#ifndef BAM_OUTFILE_HPP
#define BAM_OUTFILE_HPP

#include <string>

#include <htslib/sam.h>

#include "bam_record.hpp"

class bam_outfile {
public:
  bam_outfile(): fp(nullptr), hdr(nullptr), error_code(0) {}

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

  bam_outfile(const std::string &filename, const sam_hdr_t *const init_hdr,
              const bool bam_fmt = false) {
    if (init_hdr) {
      fp = hts_open(filename.c_str(), (bam_fmt ? "bw" : "w"));
      error_code = (fp) ? 0 : -1;
      if (!error_code) {
        hdr = sam_hdr_dup(init_hdr);
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

  bam_outfile &put_bam_rec(const bam_rec &br) {
    if (!error_code) { // already in error; do nothing
      error_code = sam_write1(fp, hdr, &br.get_rec());
      if (error_code > 0) error_code = 0;
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
