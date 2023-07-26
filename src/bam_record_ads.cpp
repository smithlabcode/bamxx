// from HTSlib
#include <htslib/sam.h>

#include "bam_record_ads.hpp"

#include <algorithm>
#include <cstdarg>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::runtime_error;
using std::string;
using std::stringstream;
using std::vector;

using std::copy;
using std::copy_backward;
using std::copy_n;
using std::fill_n;

using std::cerr;
using std::endl;




const uint16_t freverse = 16;

template<typename T> bool
high_bit_set(T x) {
  return (x >> (sizeof(T) * 8 - 1 - std::is_signed<decltype(x)>::value)) & 1;
}

/*! @hideinitializer
  @abstract  Round up to next power of two
  @discussion
  This macro will work for unsigned types up to uint64_t.

  If the next power of two does not fit in the given type, it will set
  the largest value that does.
 */
#define x_kroundup64(x)                                                        \
  ((x) > 0 ? (--(x), (x) |= (x) >> (sizeof(x) / 8),                            \
              (x) |= (x) >> (sizeof(x) / 4), (x) |= (x) >> (sizeof(x) / 2),    \
              (x) |= (x) >> (sizeof(x)), (x) |= (x) >> (sizeof(x) * 2),        \
              (x) |= (x) >> (sizeof(x) * 4), (x) += !high_bit_set(x), (x))     \
           : 0)

// From sam.c in htslib. Seems not used anywhere and not in any interface.
static int
sam_realloc_bam_data(bam1_t *b, size_t desired) {
  uint32_t new_m_data = desired;
  x_kroundup64(new_m_data);
  if (new_m_data < desired) {
    errno = ENOMEM; // (from sam.c) Not strictly true but we can't
                    // store the size
    return -1;
  }
  uint8_t *new_data = nullptr;
  if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
    new_data = static_cast<uint8_t *>(realloc(b->data, new_m_data));
  }
  else {
    new_data = static_cast<uint8_t *>(malloc(new_m_data));
    if (new_data != nullptr) {
      if (b->l_data > 0)
        copy_n(new_data, (static_cast<uint32_t>(b->l_data) < b->m_data) ? 
            b->l_data : b->m_data, b->data);
      bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
    }
  }
  if (!new_data) return -1;
  b->data = new_data;
  b->m_data = new_m_data;
  return 0;
}

static inline void
set_qseq_internal(const string &s, bam1_t *b) {
  const size_t prev_l_data = b->l_data;
  // this includes both the sequence and the qual part
  auto data_end = b->data + prev_l_data;
  const int32_t l_qseq = b->core.l_qseq;
  const int32_t delta_bytes =
    ((s.size() + 1) / 2 - (l_qseq + 1) / 2) + (s.size() - l_qseq);
  if (delta_bytes > 0) { // more memory needed
    const size_t updated_l_data = prev_l_data + delta_bytes;
    // this below handles memory policy and m_data
    const int err = sam_realloc_bam_data(b, updated_l_data);
    if (err) throw runtime_error("failed to allocate in set_seq");
    // move the data after the qseq to its new place backwards
    copy_backward(bam_get_aux(b), data_end, data_end + delta_bytes);
  }
  else if (delta_bytes < 0) // less memory needed
    copy(bam_get_aux(b), data_end, bam_get_aux(b) + delta_bytes);

  // fill in the updated query sequence
  auto dat_itr = bam_get_seq(b);
  size_t i = 0;
  for (; i + 1 < s.size(); i += 2)
    *dat_itr++ = (seq_nt16_table[static_cast<uint8_t>(s[i] << 4)]) | 
      seq_nt16_table[static_cast<uint8_t>(s[i + 1])];
  if (i < s.size()) *dat_itr++ = seq_nt16_table[static_cast<uint8_t>(s[i])] << 4;

  b->core.l_qseq = s.size(); // assign the query sequence length
  b->l_data += delta_bytes;  // update total data size
}

void
bam_rec::set_seq(const string &s) {
  set_qseq_internal(s, record);
  fill_n(get_qual(), s.size(), '\xff');
}

void
bam_rec::set_seq(const string &s, const string &q) {
  set_qseq_internal(s, record);
  copy(begin(q), end(q), get_qual());
}

void
bam_rec::set_cigar(const size_t cig_len, const uint32_t *cig) {
  const size_t prev_l_data = record->l_data;
  const auto data_end = record->data + prev_l_data;
  const int delta_bytes = (cig_len - n_cigar()) * sizeof(uint32_t);
  if (delta_bytes > 0) { // more memory needed
    const size_t updated_l_data = prev_l_data + delta_bytes;
    const int err = sam_realloc_bam_data(record, updated_l_data);
    if (err) throw runtime_error("failed to allocate in set_cigar");
    // move the data after the cigar to its new place backwards
    copy_backward(get_qseq(), data_end, data_end + delta_bytes);
  }
  else if (delta_bytes < 0) // less memory needed
    copy(get_qseq(), data_end, get_qseq() + delta_bytes);

  record->core.n_cigar = cig_len;
  record->l_data += delta_bytes;

  // fill in the updated cigar
  std::copy(cig, cig + cig_len, get_cigar());
}

bool
bam_rec::validate(const sam_hdr_t *const hdr) const {
  // ensure that tid is valid
  if (!(record->core.tid < hdr->n_targets)) return false;
  // TODO: check SAM standard to ensure that mtid is valid

  // ensure that the mapping location is consistent with the rlen
  if (!(record->core.pos < sam_hdr_tid2len(hdr, get_tid()))) return false;
  if (!(record->core.mpos < sam_hdr_tid2len(hdr, get_mtid()))) return false;

  // cigar and qseq are consistent
  if (!(static_cast<size_t>(get_l_qseq()) == qlen_from_cigar())) return false;

  // ensure that the space allocated for qual is consistent with qlen;
  // this can't be done directly without examining the possible
  // letters in the qseq portion and identifying the validity of the
  // start of the aux portion.
  if (false) return false;

  // ensure that the alphabets or types are valid for the initial parts

  return true;
}

// MN: commented below since I coult not find corresponding function.
//string
//bam_rec::tostring(const bam_header &hdr) const {
  //return tostring(hdr.h);
//}

bool
bam_rec::is_rev() const {
  return (record->core.flag & freverse) != 0;
}

bool
bam_rec::is_a_rich() const { 
  return bam_aux2A(bam_aux_get(record, "CV")) == 'A'; 
}

string
bam_rec::tostring(const sam_hdr_t *const hdr) const {
  kstring_t s = KS_INITIALIZE;
  // ADS: `sam_format1` passes to `sam_format_append` which returns -1
  // on error or the length of the kstring_t `s` argument
  // otherwise. It sets errno=EINVAL on "corrupted data" or
  // errno=ENOMEM on failure to grow s on appending.

  // ===> an error_code of 0 below should mean nothing happened, which
  // is an error for us.
  const int error_code = sam_format1(hdr, record, &s);
  if (error_code <= 0) { // ** see above
    if (s.s) {
      free(s.s);
      s.s = nullptr;
    }
    throw std::runtime_error("fail bam_rec::tostring");
  }
  const string r(s.s); // copy this before we free
  free(s.s);
  return r;
}


/* This table converts 2 bases packed in a byte to their reverse
 * complement. The input is therefore a unit8_t representing 2 bases.
 * It is assumed that the input uint8_t value is of form "xx" or "x-",
 * where 'x' a 4-bit number representing either A, C, G, T, or N and
 * '-' is 0000.  For example, the ouptut for "AG" is "CT". The format
 * "x-" is often used at the end of an odd-length sequence.  The
 * output of "A-" is "-T", and the output of "C-" is "-G", and so
 * forth. The user must handle this case separately.
 */
const uint8_t byte_revcomp_table[] = {
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  8, 136, 72, 0, 40, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 248,
  4, 132, 68, 0, 36, 0, 0, 0, 20, 0, 0, 0, 0, 0, 0, 244,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  2, 130, 66, 0, 34, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 242,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  1, 129, 65, 0, 33, 0, 0, 0, 17, 0, 0, 0, 0, 0, 0, 241,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
  0,   0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,   0,
 15, 143, 79, 0, 47, 0, 0, 0, 31, 0, 0, 0, 0, 0, 0, 255
};


static inline void
revcomp_byte_then_reverse(uint8_t *a, uint8_t *b) {
  uint8_t *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 = byte_revcomp_table[*p1];
    *p2 = byte_revcomp_table[*p2];
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
  }
  if (p1 == p2) *p1 = byte_revcomp_table[*p1];
}

void
bam_rec::revcomp() {
  const size_t l_qseq = bam_rec::get_l_qseq();
  auto seq = bam_rec::get_qseq();
  const size_t num_bytes = (l_qseq + 1) / 2;
  auto seq_end = seq + num_bytes;
  revcomp_byte_then_reverse(seq, seq_end);
  if (l_qseq % 2 == 1) { // for odd-length sequences
    for (size_t i = 0; i < num_bytes - 1; i++) {
      // swap 4-bit chunks within consecutive bytes like this:
      // (----aaaa bbbbcccc dddd....) => (aaaabbbb ccccdddd ....)
      seq[i] = (seq[i] << 4) | (seq[i + 1] >> 4);
    }
    seq[num_bytes - 1] <<= 4;
  }
}



void
bam_rec::flip() {
  if (bam_rec::is_rev())
    record->core.flag = record->core.flag & (~BAM_FREVERSE);
  else
    record->core.flag = record->core.flag | BAM_FREVERSE;

  bam_rec::revcomp();

  // ADS: don't like *(cv + 1) below, but no HTSlib function for it?
  uint8_t *cv = bam_aux_get(record, "CV");
  if (!cv) throw runtime_error("bam_aux_get failed for CV");
  *(cv + 1) = 'T';
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
bam_set1_wrapper(bam1_t *bam,
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
  bam_set1_core(bam->core, l_qname, flag, tid, pos, mapq, n_cigar,
                mtid, mpos, isize, l_seq, qname_nuls);

  const size_t data_len =
    (l_qname + qname_nuls + n_cigar*sizeof(uint32_t) + (l_seq + 1) / 2 + l_seq);

  bam->l_data = data_len;
  if (data_len + l_aux > bam->m_data) {
    const int ret = sam_realloc_bam_data(bam, data_len + l_aux);
    if (ret < 0) {
      throw runtime_error("Failed to allocate memory for BAM record");
    }
  }
  auto data_iter = bam->data;

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




