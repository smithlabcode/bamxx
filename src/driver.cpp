
#include "bam_record.hpp"

// from HTSlib
#include <htslib/sam.h>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "bam_record.hpp"
#include "uniq.hpp"


using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::cerr;
using std::sort;


size_t bam_rec_count = 0;

int main(const int argc, const char *argv[]) {

  //string inputfile(argv[1]);
  //string outfile(argv[2]);

  //bam_infile bf(inputfile);
  //if (bf.is_bam_or_sam()) {
    //cout << "bam or sam" << endl;
  //}

  //bam_header bh(bf.file->bam_header, false);
  //bam_outfile bo(outfile, bh);


  //size_t N = 5;
  //size_t count = 0;

  //bam_rec aln;
  //vector<bam_rec> bam_vec;
  //while ((bf >> aln) && count < N) {
    //bam_vec.push_back(aln); 
    //count++;
  //}
  //cout << "start sorting" << endl;
  //sort(begin(bam_vec), end(bam_vec), precedes_by_end_and_strand);
  //cout << "end sorting" << endl;

  //bool precedes = precedes_by_start(aln, aln2);
  //if (precedes) cout << "a precedes" << endl;
  //else cout << "b precedes " << endl;

  //if (equivalent_chrom_and_start(aln, aln2)) cout << "equivalent" << endl;
  //else cout << "not equivalent" << endl;

  //rd_stats stats; 
  //rd_stats stats_out; 


  //size_t count = 0;
  //while (bf >> aln && count < 100) {
    //stats.update(aln);
    ////cout << aln.tostring() << endl;
    //string aln_str = aln.tostring();
    //size_t qlen = aln.qlen_from_cigar();
    //qlen++;
    //string qname = aln.qname();
    //if (bo << aln) {
      //stats_out.update(aln);
    //}
    //count += 1;
  //}
  //size_t reads_duped = 3;
  //string statfile = "statfile.txt";

  //write_stats_output(stats, stats_out, reads_duped, statfile);

  string uniq_infile(argv[1]);
  string uniq_outfile(argv[2]);

  bool verbose = true;
  size_t n_threads = 1;
  string cmd = "some command";
  string statsfile = "stats.txt";
  string histfile = "hist.txt";
  bool bam_format = true;

  uniq(verbose, n_threads, cmd, uniq_infile, statsfile, histfile, 
      bam_format, uniq_outfile);

  //cout << bf.error_code << endl
       //<< count << endl;
  return 0;
}
