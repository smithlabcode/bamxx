
#ifndef BAM_TPOOL_HPP
#define BAM_TPOOL_HPP

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "bam_infile.hpp"
#include "bam_outfile.hpp"


class bam_tpool{
public:
  bam_tpool(const int n, const int qsize) {
    tpool.pool = hts_tpool_init(n);
    tpool.qsize = qsize;
  }

  ~bam_tpool() {
    hts_tpool_destroy(tpool.pool);
  }

  template<class T> int
  set_io(const T &bam_io_file) {
    return(hts_set_thread_pool(bam_io_file.get_fp(), &tpool));
  }

  htsThreadPool tpool;
};

#endif
