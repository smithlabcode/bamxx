
#include "bam_record.hpp"

// from HTSlib
#include <htslib/sam.h>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "bam_record.hpp"


using std::string;
using std::cout;
using std::endl;
using std::vector;




int main(const int argc, const char *argv[]) {

  string inputfile(argv[1]);
  string outfile(argv[2]);

  bam_infile bf(inputfile);
  bam_rec aln;

  bam_header bh(bf.file->bam_header);

  bam_outfile bo(outfile, bh);

  size_t count = 0;
  while (bf >> aln && count < 100) {
    cout << "main=" << aln.tostring() << endl;
    bo << aln;
    count += 1;
  }
  cout << bf.error_code << endl
       << count << endl;
}
