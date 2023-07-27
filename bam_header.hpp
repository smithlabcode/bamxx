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
