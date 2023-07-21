// from HTSlib
#include <htslib/sam.h>

#include "bam_record.hpp"

#include <cstdarg>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

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






/*  Functions needed for bam_set1_wrapper begin */
static inline void
roundup_to_power_of_2(uint32_t &x) {
  bool k_high_bit_set = (x >> (sizeof(uint32_t) * 8 - 1)) & 1;
  if (x > 0) {
    uint8_t size = sizeof(uint32_t);
    --x;
    x |= x >> (size / 4);
    x |= x >> (size / 2);
    x |= x >> (size);
    x |= x >> (size * 2);
    x |= x >> (size * 4);
    x += !k_high_bit_set;
  }
  else x = 0;
}

static int
sam_realloc_bam_data(bam_rec &b, size_t desired) {
  /* returns flag: either 0 for success or -1 for error (unable to
     allocate desired memory) */
  uint32_t new_m_data = desired;
  roundup_to_power_of_2(new_m_data);
  if (new_m_data < desired)  return -1;
  uint8_t *new_data = (uint8_t *)realloc(b.record->data, new_m_data);
  if (!new_data) return -1;
  // ADS: what would be the state of members below if -1 was returned?
  b.record->data = new_data;
  b.record->m_data = new_m_data;
  return 0;
}

static inline void
bam_set1_core(bam1_core_t &core,
              const size_t l_qname, const uint16_t flag, const int32_t tid,
              const hts_pos_t pos, const uint8_t mapq, const size_t n_cigar,
              const int32_t mtid, const hts_pos_t mpos, const hts_pos_t isize,
              const size_t l_seq, const size_t qname_nuls) {
  /* ADS: These are used in `hts_reg2bin` from `htslib/hts.h` and
     likely mean "region to bin" for indexing */
  /* MN: hts_reg2bin categorizes the size of the reference region.
     Here, we use the numbers used in htslib/cram/cram_samtools.h */
  static const int min_shift = 14;
  static const int n_lvls = 5;

  core.pos = pos;
  core.tid = tid;
  core.bin = hts_reg2bin(pos, pos + isize, min_shift, n_lvls);
  core.qual = mapq;
  core.l_extranul = qname_nuls - 1;
  core.flag = flag;
  core.l_qname = l_qname + qname_nuls;
  core.n_cigar = n_cigar;
  core.l_qseq = l_seq;
  core.mtid = mtid;
  core.mpos = mpos;
  core.isize = isize;
}

int
bam_set1_wrapper(bam_rec &bam,
                 const size_t l_qname, const char *qname,
                 const uint16_t flag, const int32_t tid,
                 const hts_pos_t pos, const uint8_t mapq,
                 const size_t n_cigar, const uint32_t *cigar,
                 const int32_t mtid, const hts_pos_t mpos,
                 const hts_pos_t isize, const size_t l_seq,
                 const size_t l_aux) {
  /* This is based on how assignment is done in the `bam_set1`
     function defined in `sam.c` from htslib */

  /*
   * This modification assigns variables of bam1_t struct but not the sequence.
   *
   * Many checks have been removed because they are checked in code
   * that calls this function, mostly because they already come from a
   * valid `bam1_t` struct and so the values have been individually
   * validated.
   *
   * Assumptions:
   * cigar has been computed and is in the right format
   * rlen = isize
   * qlen = l_seq
   * l_qname <= 254
   * HTS_POS_MAX - rlen > pos
   * Where HTS_POS_MAX = ((((int64_t)INT_MAX)<<32)|INT_MAX) is the highest
   * supported position.
   *
   * Number of bytes needed for the data is smaller than INT32_MAX
   *
   * qual = NULL, because we do not keep the quality scores through
   * formatting the reads.
   */

  // `qname_nuls` below is the number of '\0' to use to pad the qname
  // so that the cigar has 4-byte alignment.
  const size_t qname_nuls = 4 - l_qname % 4;
  bam_set1_core(bam.record->core, l_qname, flag, tid, pos, mapq, n_cigar,
                mtid, mpos, isize, l_seq, qname_nuls);

  const size_t data_len =
    (l_qname + qname_nuls + n_cigar*sizeof(uint32_t) + (l_seq + 1) / 2 + l_seq);

  bam.record->l_data = data_len;
  if (data_len + l_aux > bam.record->m_data) {
    const int ret = sam_realloc_bam_data(bam, data_len + l_aux);
    if (ret < 0) {
      throw std::runtime_error("Failed to allocate memory for BAM record");
    }
  }
  auto data_iter = bam.record->data;

  std::copy_n(qname, l_qname, data_iter);
  std::fill_n(data_iter + l_qname, qname_nuls, '\0');
  data_iter += l_qname + qname_nuls;

  // ADS: reinterpret here because we know the cigar is originally an
  // array of uint32_t and has been aligned for efficiency
  std::copy_n(cigar, n_cigar, reinterpret_cast<uint32_t *>(data_iter));
  data_iter += n_cigar * sizeof(uint32_t);

  // skipping sequece assignment
  data_iter += (l_seq + 1) / 2;

  std::fill(data_iter, data_iter + l_seq, '\xff');

  return static_cast<int>(data_len);
}


/*  Functions needed for bam_set1_wrapper end */


