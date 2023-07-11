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


inline
void write_cigar(std::stringstream &ss, const bam_rec& br) {
  uint32_t* cigar = bam_get_cigar(br.get());
  for (size_t i = 0; i < br.n_cigar(); i++) {
    ss << bam_cigar_oplen(*cigar) << bam_cigar_opchr(*cigar);
    cigar++;
  }
  ss << '\t';
}

inline 
void write_seq(std::stringstream &ss, const bam_rec& br) {
  uint8_t* seq = bam_get_seq(br.get());
  for (auto i = 0; i < br.l_qseq(); i++) {
    ss << seq_nt16_str[bam_seqi(seq, i)];
  }
  ss << '\t';
}

inline 
void write_qual(std::stringstream &ss, const bam_rec& br) {
  if (br.qual() == 0xff) ss << "*";
  else {
    uint8_t* qual_seq = bam_get_qual(br.get());
    for (auto i = 0; i < br.l_qseq(); i++) {
      ss << static_cast<char>(qual_seq[i]);
    }
  }
}

inline 
void write_aux(std::stringstream &ss, const bam_rec& br) {
  kstring_t str = {0,0, NULL};
  uint8_t* aux = bam_get_aux(br.get()); 
  const uint8_t* end = br.data() + br.l_data();
  while (end - aux >= 4) {
    kputc_('\t', &str);
    aux = const_cast<uint8_t*>
            (sam_format_aux1(aux, aux[2], aux+3, end, &str)); 
  }
  kputsn("", 1, &str); 
  ss << str.s;
  free(str.s);
}

inline 
void write_mtid(std::stringstream &ss, const bam_rec& br) {
  if (br.mtid() == -1) ss << "*\t";
  else ss << br.mtid() << "\t";
}

string 
bam_rec::tostring() const {
  std::stringstream ss;

  ss << bam_get_qname(record) << "\t";
  ss << flag() << "\t";
  ss << tid() << "\t";
  ss << pos()+1 << "\t";
  ss << static_cast<int>(qual()) << "\t";
  write_cigar(ss, *this);
  write_mtid(ss, *this);
  ss << mpos()+1 << "\t";
  ss << isize() << "\t";
  write_seq(ss, *this);
  write_qual(ss, *this);
  write_aux(ss,*this);

  return ss.str();
}

