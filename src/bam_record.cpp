// from HTSlib
#include <htslib/sam.h>


#include "bam_record.hpp"


#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::cout;
using std::endl;
using std::vector;



string 
bam_rec::tostring() const {
  std::stringstream ss;
  ss << bam_get_qname(record) << "\t";
  ss << this->flag() << "\t";
  ss << this->tid() << "\t";
  ss << this->pos()+1 << "\t";
  ss << static_cast<int>(this->qual()) << "\t";

  uint32_t* cigar = bam_get_cigar(record);
  for (size_t i = 0; i < this->n_cigar(); i++) {
    ss << bam_cigar_oplen(*cigar) << bam_cigar_opchr(*cigar);
    cigar++;
  }
  ss << "\t";

  if (this->mtid() == -1) ss << "*\t";
  else ss << this->mtid() << "\t";

  ss << this->mpos()+1 << "\t";
  ss << this->isize() << "\t";

  uint8_t* seq = bam_get_seq(record);
  for (size_t i = 0; i < static_cast<size_t>(this->l_qseq()); i++) {
    ss << seq_nt16_str[bam_seqi(seq, i)];
  }
  ss << "\t";
  
  if (this->qual() == 0xff) ss << "*";
  else {
    uint8_t* qual = bam_get_qual(record);
    for (size_t i = 0; i < static_cast<size_t>(this->l_qseq()); i++) {
      ss << static_cast<char>(qual[i]);
    }
  }
  
  kstring_t str = {0,0, NULL};
  uint8_t* aux = bam_get_aux(record); 
  uint8_t* end = record->data + record->l_data;
  while (end - aux >= 4) {
    kputc_('\t', &str);
    aux = const_cast<uint8_t*>(sam_format_aux1(aux, aux[2], aux+3, end, &str)); 
  }
  kputsn("", 1, &str); 
  ss << str.s;
  free(str.s);

  return ss.str();
}

