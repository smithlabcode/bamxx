// from HTSlib
#include <htslib/sam.h>

#include "bam_record.hpp"

#include <cstdarg>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;


const uint8_t cigar_shift = 4;
const uint32_t cigar_mask = 0xf;
const uint8_t sequence_data_enum = 1;
const uint8_t sam_enum = 3;
const uint8_t bam_enum = 4;
const uint16_t freverse = 16;



static inline int
str_resize(kstring_t *str, size_t size) {
  if (str->m < size) {
    // size = (size > (SIZE_MAX>>2)) ? size : size + (size >> 1);
    size = (size > (SIZE_MAX / 4)) ? size : size * 3;
    char *tmp = static_cast<char *>(realloc(str->s, size));
    if (!tmp) return -1;
    str->s = tmp;
    str->m = size;
  }
  return 0;
}

static inline int
putsn(const char *p, size_t l, kstring_t *s) {
  size_t new_sz = s->l + l + 2;
  if (new_sz <= s->l || str_resize(s, new_sz) < 0) return EOF;
  memcpy(s->s + s->l, p, l);
  s->l += l;
  s->s[s->l] = 0;
  return l;
}

static inline int
putc(int c, kstring_t *s) {
  if (str_resize(s, s->l + 1) < 0) return EOF;
  s->s[s->l++] = c;
  return 1;
}

inline size_t
cigar_oplen(const uint32_t c) {
  return c >> cigar_shift;
}

inline uint8_t
cigar_op(const uint32_t c) {
  return c & cigar_mask;
}

static const char cigar_str[] = "MIDNSHP=XB??????";

inline char
cigar_opchr(const uint32_t c) {
  return cigar_str[cigar_op(c)];
}


void
bam_rec::write_cigar(stringstream &ss) const {
  const uint32_t *cigar = bam_rec::get_cigar();
  for (size_t i = 0; i < record->core.n_cigar; i++) {
    ss << cigar_oplen(*cigar) << cigar_opchr(*cigar);
    cigar++;
  }
  ss << '\t';
}



static const char *const byte2str[] = {
  "",   "A",  "C",  "NN", "G",  "NN", "NN", "NN", "T",  "NN", "NN", "NN", "NN",
  "NN", "NN", "N",  "A",  "AA", "AC", "NN", "AG", "NN", "NN", "NN", "AT", "NN",
  "NN", "NN", "NN", "NN", "NN", "AN", "C",  "CA", "CC", "NN", "CG", "NN", "NN",
  "NN", "CT", "NN", "NN", "NN", "NN", "NN", "NN", "CN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "G",
  "GA", "GC", "NN", "GG", "NN", "NN", "NN", "GT", "NN", "NN", "NN", "NN", "NN",
  "NN", "GN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "T",  "TA",
  "TC", "NN", "TG", "NN", "NN", "NN", "TT", "NN", "NN", "NN", "NN", "NN", "NN",
  "TN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "N",  "NA", "NC", "NN", "NG", "NN", "NN",
  "NN", "NT", "NN", "NN", "NN", "NN", "NN", "NN", "NN"
};


void
bam_rec::write_seq(stringstream &ss) const {
  const uint8_t *seq = bam_rec::get_seq();
  for (auto i = 0; i < (record->core.l_qseq + 1) / 2; i++) { 
    ss << byte2str[seq[i]]; 
  }
  ss << '\t';
}

void
bam_rec::write_qual(stringstream &ss) const {
  if (record->core.qual == 0xff)
    ss << "*";
  else {
    const uint8_t *qual_seq = bam_rec::get_qual();
    for (auto i = 0; i < record->core.l_qseq; i++) {
      ss << static_cast<char>(qual_seq[i]);
    }
  }
}

void
bam_rec::write_aux(std::stringstream &ss) const {
  kstring_t str = {0, 0, NULL};
  const uint8_t *aux = bam_rec::get_aux();
  const uint8_t *end = record->data + record->l_data;
  while (end - aux >= 4) {
    putc('\t', &str);
    aux = static_cast<const uint8_t *> // Not sure how else to do it.
      (sam_format_aux1(aux, aux[2], aux + 3, end, &str));
  }
  putsn("", 0, &str);
  ss << str.s;
  free(str.s);
}

inline void
bam_rec::write_mtid(stringstream &ss) const {
  if (record->core.mtid == -1)
    ss << "*\t";
  else
    ss << record->core.mtid << "\t";
}

string
bam_rec::tostring() const {
  std::stringstream ss;
  bam1_core_t *core = &record->core;

  ss << bam_rec::qname() << "\t";
  ss << core->flag << "\t";
  ss << core->tid << "\t";
  ss << core->pos + 1 << "\t";
  ss << static_cast<int>(core->pos) << "\t";
  bam_rec::write_cigar(ss);
  bam_rec::write_mtid(ss);
  ss << core->mpos + 1 << "\t";
  ss << core->isize << "\t";
  bam_rec::write_seq(ss);
  bam_rec::write_qual(ss);
  bam_rec::write_aux(ss);

  return ss.str();
}


bool
bam_rec::is_rev() const {
  return (record->core.flag & freverse) != 0;
}

// MN: Not working as I want it to.
//     Will work in the future
int
bam_header::add_line(const string &type, const vector<string> &args) {
  //va_list args;
  //va_start(args, type);
  //int err_num = sam_hdr_add_line(header, type.c_str(), args);
  //va_end(args);
  //return err_num;
  return 0;
}





/*                        bam_infile                    */ 

bool
bam_infile::is_bam_or_sam() {
  return fmt->category == sequence_data_enum &&
         (fmt->format == bam_enum || fmt->format == sam_enum);
}
