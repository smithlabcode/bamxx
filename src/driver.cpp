
#include "bam_record.hpp"

// from HTSlib
#include <htslib/sam.h>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "bam_record.hpp"
#include "uniq.hpp"


using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::cerr;




int main(const int argc, const char *argv[]) {

  string inputfile(argv[1]);
  string outfile(argv[2]);

  bam_infile bf(inputfile);
  if (bf.is_bam_or_sam()) {
    cout << "bam or sam" << endl;
  }

  bam_rec aln, aln2;

  

  bam_header bh(bf.file->bam_header);

  bam_outfile bo(outfile, bh);

  bf >> aln;
  bf >> aln2;

  bool precedes = precedes_by_start(aln, aln2);
  if (precedes) cout << "a precedes" << endl;
  else cout << "b precedes " << endl;

  if (equivalent_chrom_and_start(aln, aln2)) cout << "equivalent" << endl;
  else cout << "not equivalent" << endl;

  rd_stats stats; 
  rd_stats stats_out; 


  size_t count = 0;
  while (bf >> aln && count < 100) {
    stats.update(aln);
    //cout << aln.tostring() << endl;
    string aln_str = aln.tostring();
    size_t qlen = aln.qlen_from_cigar();
    qlen++;
    string qname = aln.qname();
    if (bo << aln) {
      stats_out.update(aln);
    }
    count += 1;
  }
  size_t reads_duped = 3;
  string statfile = "statfile.txt";

  write_stats_output(stats, stats_out, reads_duped, statfile);


  bool verbose = true;
  size_t n_threads = 1;
  string cmd = "some command";
  string uniq_infile = "test_data/reads.fmt.bam";
  string statsfile = "stats.txt";
  string histfile = "hist.txt";
  bool bam_format = true;
  string uniq_outfile = "uniq.out.sam";

  //uniq(verbose, n_threads, cmd, uniq_infile, statsfile, histfile, 
      //bam_format, uniq_outfile);





  cout << bf.error_code << endl
       << count << endl;
}
