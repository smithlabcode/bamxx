
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


  size_t count = 0;
  while (bf >> aln && count < 100) {
    cout << aln.tostring() << endl;
    size_t qlen = aln.qlen_from_cigar();
    qlen++;
    string qname = aln.qname();
    bo << aln;
    count += 1;
  }
  cout << bf.error_code << endl
       << count << endl;
}
