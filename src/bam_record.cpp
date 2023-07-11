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
  ss << flag() << "\t";
  ss << tid() << "\t";
  ss << pos()+1 << "\t";
  ss << static_cast<int>(qual()) << "\t";

  uint32_t* cigar = bam_get_cigar(record);
  for (size_t i = 0; i < n_cigar(); i++) {
    ss << bam_cigar_oplen(*cigar) << bam_cigar_opchr(*cigar);
    cigar++;
  }
  ss << "\t";

  if (mtid() == -1) ss << "*\t";
  else ss << mtid() << "\t";

  ss << mpos()+1 << "\t";
  ss << isize() << "\t";

  uint8_t* seq = bam_get_seq(record);
  for (auto i = 0; i < l_qseq(); i++) {
    ss << seq_nt16_str[bam_seqi(seq, i)];
  }
  ss << "\t";
  
  if (qual() == 0xff) ss << "*";
  else {
    uint8_t* qual = bam_get_qual(record);
    for (auto i = 0; i < l_qseq(); i++) {
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

