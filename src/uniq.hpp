
#ifndef UNIQ_HPP
#define UNIQ_HPP

#include <string>
#include <htslib/sam.h>

#include <iostream>
#include <cassert>
#include <vector>

#include "bam_record.hpp"


inline bool
precedes_by_start(const bam_rec &a, const bam_rec &b) {
  return a.tid() == b.tid() && a.pos() < b.pos();
}


inline bool
equivalent_chrom_and_start(const bam_rec &a, const bam_rec &b) {
  return a.pos() == b.pos() && a.tid() == b.tid();
}

struct rd_stats { // keep track of good bases/reads in and out
  size_t bases;
  size_t reads;
  rd_stats() : bases(0), reads(0) {}
  void update(const bam_rec &br) { 
    //MN: it might be faster to use br.l_qseq() if we know
    //    br is well defined.
    bases += br.qlen_from_cigar(); 
    ++reads;
  }
};


void
write_stats_output(const rd_stats &rs_in, const rd_stats &rs_out,
                   const size_t reads_duped, const std::string &statfile);


// todo: test in drive once vector is ready
void
write_hist_output(const std::vector<size_t> &hist, 
    const std::string &histfile);



#endif
