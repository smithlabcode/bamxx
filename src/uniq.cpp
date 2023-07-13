
// from HTSlib
#include <htslib/sam.h>


#include "uniq.hpp"


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::runtime_error;
using std::ofstream;


void
write_stats_output(const rd_stats &rs_in, const rd_stats &rs_out,
                   const size_t reads_duped, const string &statfile) {
  if (!statfile.empty()) {
    const size_t reads_removed = rs_in.reads - rs_out.reads;
    const double non_dup_frac =
      (rs_out.reads - reads_duped)/static_cast<double>(rs_in.reads);
    const double dup_rate =
      (reads_removed + reads_duped)/static_cast<double>(reads_duped);
    ofstream out_stat(statfile);
    if (!out_stat) throw runtime_error("bad stats output file");
    out_stat << "total_reads: " << rs_in.reads << endl
             << "total_bases: " << rs_in.bases << endl
             << "unique_reads: " << rs_out.reads << endl
             << "unique_read_bases: " << rs_out.bases << endl
             << "non_duplicate_fraction: " << non_dup_frac << endl
             << "duplicate_reads: " << reads_duped << endl
             << "reads_removed: " << reads_removed << endl
             << "duplication_rate: " << dup_rate << endl;
  }
}

void
write_hist_output(const vector<size_t> &hist, const string &histfile) {
  if (!histfile.empty()) {
    ofstream out_hist(histfile);
    if (!out_hist) throw runtime_error("bad hist output file");
    for (size_t i = 0; i < hist.size(); ++i)
      if (hist[i] > 0)
        out_hist << i << '\t' << hist[i] << '\n';
  }
}



/* The "inner" buffer corresponds to all reads sharing chrom, start,
 * end and strand, and is a contiguous subset of the "outer" buffer
 * that shares the same end and strand.
 */
//static void
//process_inner_buffer(const vector<bam_rec>::const_iterator it,
                     //const vector<bam_rec>::const_iterator jt,
                     //bam_header &hdr, bam_outfile &out,
                     //rd_stats &rs_out,
                     //size_t &reads_duped,
                     //vector<size_t> &hist) {
  //const size_t n_reads = std::distance(it, jt);
  //const size_t selected = rand() % n_reads;
  //if (sam_write1(out, hdr, *(it + selected)) < 0)
    //throw runtime_error("failed writing bam record");
  //if (hist.size() <= n_reads)
    //hist.resize(n_reads + 1);
  //hist[n_reads]++;
  //rs_out.update(*(it + selected));
  //reads_duped += (n_reads > 1);
//}
