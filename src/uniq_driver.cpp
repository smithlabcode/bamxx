
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

  string inputfile(argv[1]);
  string outfile(argv[2]);

  bam_infile bf(inputfile);

  size_t N = 5;
  size_t count = 0;

  bam_rec aln;
  vector<bam_rec> bam_vec;
  while ((bf >> aln) && count < N) {
    bam_vec.push_back(aln); 
    count++;
  }
  cout << "start sorting" << endl;
  sort(begin(bam_vec), end(bam_vec), precedes_by_end_and_strand);
  cout << "end sorting" << endl;
  return 0;
}
