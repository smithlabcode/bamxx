
// from HTSlib
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "uniq.hpp"


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::runtime_error;
using std::ofstream;
using std::sort;
using std::to_string;


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
void
process_inner_buffer(const vector<bam_rec>::const_iterator it,
                     const vector<bam_rec>::const_iterator jt,
                     bam_header &hdr, bam_outfile &bo,
                     rd_stats &rs_out,
                     size_t &reads_duped,
                     vector<size_t> &hist) {
  const size_t n_reads = std::distance(it, jt);
  const size_t selected = rand() % n_reads;
  
  if (!bo.put_bam_rec(*(it + selected)))
    throw runtime_error("failed writing bam record");
  if (hist.size() <= n_reads)
    hist.resize(n_reads + 1);
  hist[n_reads]++;
  rs_out.update(*(it + selected));
  reads_duped += (n_reads > 1);
}

/* The buffer corresponds to reads sharing the same mapping chromosome
   and start position. These are gathered and then processed together. */
void
process_buffer(rd_stats &rs_out, size_t &reads_duped, 
               vector<size_t> &hist, vector<bam_rec> &buffer, 
               bam_header &hdr, bam_outfile &out) {
  sort(begin(buffer), end(buffer), precedes_by_end_and_strand);
  auto it(begin(buffer));
  auto jt = it + 1;
  for (; jt != end(buffer); ++jt)
    if (!equivalent_end_and_strand(*it, *jt)) {
      process_inner_buffer(it, jt, hdr, out, rs_out, reads_duped, hist);
      it = jt;
    }
  process_inner_buffer(it, jt, hdr, out, rs_out, reads_duped, hist);
  buffer.clear();
}



void
uniq(const bool VERBOSE, const size_t n_threads,
     const string &cmd, const string &infile,
     const string &statfile, const string &histfile,
     const bool bam_format, const string &outfile) {

  bam_infile hts(infile);
  if (!hts || errno)
    throw runtime_error("bad htslib file: " + infile);

  bam_tpool thread_pool(n_threads, 0);
  if (thread_pool.set(hts) < 0)
    throw runtime_error("error setting threads");

  if (!hts.is_bam_or_sam())
    throw runtime_error("bad file format: " + infile);

  //if (!hdr)
    //throw runtime_error("failed to read header: " + infile);
  bam_header hdr(hts.file->bam_header, true); // shallow copy

  //bam_outfile out = hts_open(outfile.c_str(), bam_format ? "wb" : "w");
  bam_header hdr_out(hdr, false); // deep copy

  //MN: need to replace with wrapper function  
  //    need to replace VERSION_PLACEHOLDER with VERSION
  if (sam_hdr_add_line(hdr_out.header, "PG", "ID",
                       "DNMTOOLS", "VN", "VERSION_PLACEHOLDER", 
                       "CL", cmd.c_str(), NULL))
    throw runtime_error("failed to format header");

  // MN: Currently only outputs to sam file
  bam_outfile out(outfile, hdr_out);
  if (out.error_code) {
    throw runtime_error("failed to open out file");
  }

  if (thread_pool.set(out) < 0)
    throw runtime_error("error setting threads");

  // values to tabulate stats; no real cost
  rd_stats rs_in, rs_out;
  size_t reads_duped = 0;
  vector<size_t> hist;

  bam_rec aln; 
  if (!(hts >> aln))
    throw runtime_error("failed parsing read from input file");
  rs_in.update(aln); // update for the input we just got

  vector<bam_rec> buffer; // select output from this buffer
  buffer.push_back(aln);

  vector<bool> chroms_seen(hdr.n_targets(), false);
  int32_t cur_chrom = aln.tid();

  while (hts >> aln) {
    rs_in.update(aln);

    // below works because buffer reset at every new chrom
    if (precedes_by_start(aln, buffer[0]))
      throw runtime_error("not sorted: " + buffer[0].qname() + 
          " " + aln.qname());

    const int32_t chrom = aln.tid();
    if (chrom != cur_chrom) {
      if (chroms_seen[chrom]) throw runtime_error("input not sorted");
      chroms_seen[chrom] = true;
      cur_chrom = chrom;
    }

    if (!equivalent_chrom_and_start(buffer[0], aln))
      process_buffer(rs_out, reads_duped, hist, buffer, hdr, out);
    buffer.push_back(aln);
  }
  process_buffer(rs_out, reads_duped, hist, buffer, hdr, out);

  // write any additional output requested
  write_stats_output(rs_in, rs_out, reads_duped, statfile);
  write_hist_output(hist, histfile);
}
