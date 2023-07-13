

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
using std::cerr; 



void test_bam1_t(const string& inputfile, const string &outfile) {
  htsFile* bf = hts_open(inputfile.c_str(), "r");
  sam_hdr_t *hdr = bam_hdr_read(bf->fp.bgzf);
  bf->bam_header = hdr;
  
  htsFile* bo = hts_open(outfile.c_str(), "w");
  bo->bam_header = hdr;  
  bo->bam_header->ref_count++;
  int errnum = sam_hdr_write(bo, bo->bam_header);
  if (errnum) {
    cerr << "Failed to write header" << endl;
    std::exit(1);
  }

  bam1_t *aln = bam_init1();
  
  int err_code = 0;

  while((err_code = sam_read1(bf, hdr, aln)) >= 0) {
    aln->core.pos += 123;
    if(sam_write1(bo, bo->bam_header, aln)) {
      cerr << "Failed to write read" << endl;
      std::exit(1);
    }
  }

  bam_destroy1(aln);
  hts_close(bo);
  hts_close(bf);

}

void test_bam_rec(const string& inputfile, const string &outfile) {
  bam_infile bf(inputfile);
  bam_header bh(bf.file->bam_header);
  bam_outfile bo(outfile, bh);

  bam_rec aln;
  
  while((bf >> aln)) {
    aln.pos() += 123;
    bo << aln;
  }
}

int main(const int argc, const char *argv[]) {

  string inputfile(argv[1]);
  string outfile(argv[2]);
  string mode(argv[3]);

  if (mode == "bam1_t") {
    test_bam1_t(inputfile, outfile);
  }
  else {
    test_bam_rec(inputfile, outfile);
  }

}
