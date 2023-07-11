#include "bam_record.hpp"
// from HTSlib
#include <htslib/sam.h>

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::cout;
using std::endl;
using std::vector;



string bam_rec::tostring() const {
  std::stringstream ss;
  ss << bam_get_qname(record) << "\t";
  ss << this->flag() << "\t";
  ss << "chr" << this->tid()+1 << "\t";
  ss << this->pos()+1 << "\t";
  ss << static_cast<int>(this->qual()) << "\t";

  auto* cigar = bam_get_cigar(record);
  for (size_t i = 0; i < this->n_cigar(); i++) {
    ss << bam_cigar_oplen(*cigar) << bam_cigar_opchr(*cigar);
    cigar++;
  }
  ss << "\t";

  if (this->mtid() == -1) ss << "*\t";
  else ss << "chr" << this->mtid()+1 << "\t";

  ss << this->mpos()+1 << "\t";
  ss << this->isize() << "\t";

  auto* seq = bam_get_seq(record);
  for (size_t i = 0; i < static_cast<size_t>(this->l_qseq()); i++) {
    ss << seq_nt16_str[bam_seqi(seq, i)];
  }
  ss << "\t";
  
  if (this->qual() == 0xff) ss << "*";
  else {
    auto* qual = bam_get_qual(record);
    for (size_t i = 0; i < static_cast<size_t>(this->l_qseq()); i++) {
      ss << static_cast<char>(qual[i]);
    }
  }
  
  kstring_t str = {0,0, NULL};
  auto* aux = bam_get_aux(record); 
  auto* end = record->data + record->l_data;
  while (end - aux >= 4) {
    kputc_('\t', &str);
    aux = const_cast<uint8_t*>(sam_format_aux1(aux, aux[2], aux+3, end, &str)); 
  }
  kputsn("", 1, &str); 
  ss << str.s;
  free(str.s);

  return ss.str();
}

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