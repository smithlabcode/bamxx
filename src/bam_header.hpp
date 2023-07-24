#ifndef BAM_HEADER_HPP
#define BAM_HEADER_HPP

#include <string>
#include <stdexcept>

#include <htslib/sam.h>

struct bam_header {
  bam_header(const sam_hdr_t *const other) { h = sam_hdr_dup(other); }
  ~bam_header() { if (h != nullptr) { bam_hdr_destroy(h); h = nullptr; } }
  operator bool() const {return h != nullptr;}
  void
  sam_hdr_add_program(const std::string identifier,
                      const std::string program,
                      const std::string package,
                      const std::string version,
                      const std::string cmdline) {
    const int ret = sam_hdr_add_line(h, "PG",
                                     "ID", identifier.c_str(),
                                     "PN", program.c_str(),
                                     "PP", package.c_str(),
                                     "VN", version.c_str(),
                                     "CL", cmdline.c_str(), nullptr);
    if (ret != 0) throw std::runtime_error("fail sam_hdr_add_program");
  }
  sam_hdr_t *h{};
};

#endif
