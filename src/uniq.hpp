
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
precedes_by_end_and_strand(const bam_rec &a, const bam_rec &b) {
  const size_t end_a = a.endpos(), end_b = b.endpos();
  return (end_a < end_b || (end_a == end_b && a.is_rev() < b.is_rev()));
}

inline bool
equivalent_chrom_and_start(const bam_rec &a, const bam_rec &b) {
  return a.pos() == b.pos() && a.tid() == b.tid();
}

inline bool
equivalent_end_and_strand(const bam_rec &a, const bam_rec &b) {
  return a.endpos() == b.endpos() && a.is_rev() == b.is_rev();
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


void
process_inner_buffer(const std::vector<bam_rec>::const_iterator it,
                     const std::vector<bam_rec>::const_iterator jt,
                     bam_header &hdr, bam_outfile &bo,
                     rd_stats &rs_out,
                     size_t &reads_duped,
                     std::vector<size_t> &hist);

void
process_buffer(rd_stats &rs_out, size_t &reads_duped, 
               std::vector<size_t> &hist, std::vector<bam_rec> &buffer, 
               bam_header &hdr, bam_outfile &out);


void
uniq(const bool VERBOSE, const size_t n_threads,
     const std::string &cmd, const std::string &infile,
     const std::string &statfile, const std::string &histfile,
     const bool bam_format, const std::string &outfile);



#endif
