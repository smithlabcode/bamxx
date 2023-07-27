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

#ifndef BAM_TPOOL_HPP
#define BAM_TPOOL_HPP

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

class bam_tpool{
public:
  bam_tpool(const size_t n_threads) : tpool{hts_tpool_init(n_threads), 0} {}
  ~bam_tpool() { hts_tpool_destroy(tpool.pool); }
  template<class T> void set_io(const T &bam_file) {
    const int ret = hts_set_thread_pool(bam_file.get_fp(), &tpool);
    if (ret < 0) throw std::runtime_error("failed to set thread pool");
  }
  htsThreadPool tpool{};
};

#endif
