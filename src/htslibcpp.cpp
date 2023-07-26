
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>
#include <sstream>
#include <cstring>
#include <cstdint> // for [u]int[0-9]+_t
#include <algorithm>
#include <cassert>

#include "dnmt_error.hpp"

#include <htslib/sam.h>
#include <htslib/thread_pool.h>


// from smithlab
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

// from htslibcpp
#include "htslibcpp.hpp"


using std::string;
using std::vector;
using std::cerr;
using std::endl;





//[> Functions used in format-reads.cpp <]


bool
check_input_file(const string &infile) {
  bam_infile hts(infile);
  if (!hts || errno) throw dnmt_error("error opening: " + infile);
  const htsFormat *fmt = hts_get_format(hts.get_fp());
  if (fmt->category != sequence_data)
    throw dnmt_error("not sequence data: " + infile);
  if (fmt->format != bam && fmt->format != sam)
    throw dnmt_error("not SAM/BAM format: " + infile);

  return true;
}


vector<string>
load_read_names(const string &inputfile, const size_t n_reads) {
  bam_infile hts(inputfile);
  if (!hts) throw dnmt_error("failed to open file: " + inputfile);

  bam_rec aln;
  vector<string> names;
  size_t count = 0;

  while ((get_bam(hts, aln)) && count++ < n_reads)
    names.push_back(aln.qname());

  return names;
}

void
standardize_format(const string &input_format, bam_rec &aln) {
  int err_code = 0;

  if (input_format == "abismal" || input_format == "walt") return;

  if (input_format == "bsmap") {
    // A/T rich; get the ZS tag value
    auto zs_tag = bam_aux_get(aln.get_ptr(), "ZS");
    if (!zs_tag) throw dnmt_error("bam_aux_get for ZS (invalid bsmap)");
    // ADS: test for errors on the line below
    const uint8_t cv = string(bam_aux2Z(zs_tag))[1] == '-' ? 'A' : 'T';
    // get the "mismatches" tag
    auto nm_tag = bam_aux_get(aln.get_ptr(), "NM");
    if (!nm_tag) throw dnmt_error("bam_aux_get for NM (invalid bsmap)");
    const int64_t nm = bam_aux2i(nm_tag);

    // del aux (no data resize)
    aln.get_ptr()->l_data = bam_get_aux(aln.get_ptr()) - aln.get_ptr()->data; 

    //[> add the tags we want <]
    // "udpate" for "int" because it determines the right size; even
    // though we just deleted all tags, it will add it back here.
    err_code = bam_aux_update_int(aln.get_ptr(), "NM", nm);
    if (err_code < 0) throw dnmt_error(err_code, "bam_aux_update_int");
    // "append" for "char" because there is no corresponding update
    err_code = bam_aux_append(aln.get_ptr(), "CV", 'A', 1, &cv);
    if (err_code < 0) throw dnmt_error(err_code, "bam_aux_append");

    if (aln.is_rev())
      aln.revcomp(); // reverse complement if needed
  }
  if (input_format == "bismark") {
    // ADS: Previously we modified the read names at the first
    // underscore. Even if the names are still that way, it should no
    // longer be needed since we compare names up to a learned suffix.

    // A/T rich; get the XR tag value
    auto xr_tag = bam_aux_get(aln.get_ptr(), "XR");
    if (!xr_tag) throw dnmt_error("bam_aux_get for XR (invalid bismark)");
    const uint8_t cv = string(bam_aux2Z(xr_tag)) == "GA" ? 'A' : 'T';
    // get the "mismatches" tag
    auto nm_tag = bam_aux_get(aln.get_ptr(), "NM");
    if (!nm_tag) throw dnmt_error("bam_aux_get for NM (invalid bismark)");
    const int64_t nm = bam_aux2i(nm_tag);

    // del aux (no data resize)
    aln.get_ptr()->l_data = bam_get_aux(aln.get_ptr()) - aln.get_ptr()->data; 

    //[> add the tags we want <]
    // "udpate" for "int" because it determines the right size; even
    // though we just deleted all tags, it will add it back here.
    err_code = bam_aux_update_int(aln.get_ptr(), "NM", nm);
    if (err_code < 0) throw dnmt_error(err_code, "bam_aux_update_int");
    // "append" for "char" because there is no corresponding update
    err_code = bam_aux_append(aln.get_ptr(), "CV", 'A', 1, &cv);
    if (err_code < 0) throw dnmt_error(err_code, "bam_aux_append");

    if (bam_is_rev(aln.get_ptr()))
      aln.revcomp(); // reverse complement if needed
  }

  // Be sure this doesn't depend on mapper! Removes the "qual" part of
  // the data in a bam1_t struct but does not change its uncompressed
  // size.
  const auto qs = bam_get_qual(aln.get_ptr());
  std::fill(qs, qs + aln.get_l_qseq(), '\xff'); // deletes qseq
}

bool
same_name(const bam_rec &a, const bam_rec &b, const size_t suff_len) {
  // "+ 1" below: extranul counts *extras*; we don't want *any* nulls
  const uint16_t a_l = a.get_m_qname();
  const uint16_t b_l = b.get_m_qname();
  if (a_l != b_l) return false;
  assert(a_l > suff_len);
  return !std::strncmp(a.qname_ptr(), b.qname_ptr(), a_l - suff_len);
}

bool
are_mates(const bam_rec &one, const bam_rec &two) {
  return one.get_mtid() == two.get_tid() &&
         one.get_mpos() == two.get_mpos() &&
         one.is_rev() != two.is_rev();
}

static inline bool
eats_ref(const uint32_t c) { return bam_cigar_type(bam_cigar_op(c)) & 2; }

static inline bool
eats_query(const uint32_t c) { return bam_cigar_type(bam_cigar_op(c)) & 1; }

//MN: currenctly not used
static inline void
reverse(char *a, char *b) {
  char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
    assert(valid_base(*p1) && valid_base(*p2));
  }
}

// return value is the number of cigar ops that are fully consumed in
// order to read n_ref, while "partial_oplen" is the number of bases
// that could be taken from the next operation, which might be merged
// with the other read.

static uint32_t
get_full_and_partial_ops(const uint32_t *cig_in, const uint32_t in_ops,
                         const uint32_t n_ref_full, uint32_t *partial_oplen) {
  // assume: n_ops <= size(cig_in) <= size(cig_out)
  size_t rlen = 0;
  uint32_t i = 0;
  for (i = 0; i < in_ops; ++i) {
    if (eats_ref(cig_in[i])) {
      if (rlen + bam_cigar_oplen(cig_in[i]) > n_ref_full)
        break;
      rlen += bam_cigar_oplen(cig_in[i]);
    }
  }
  *partial_oplen = n_ref_full - rlen;
  return i;
}

static inline uint32_t
to_insertion(const uint32_t x) {
  return (x & ~BAM_CIGAR_MASK) | BAM_CINS;
}

static inline uint32_t
to_softclip(const uint32_t x) {
  return (x & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
}

static void
fix_internal_softclip(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 3) return;
  // find first non-softclip
  auto c_beg = cigar;
  auto c_end = cigar + n_cigar;

  while (!eats_ref(*c_beg) && ++c_beg != c_end);
  if (c_beg == c_end) throw dnmt_error("cigar eats no ref");

  while (!eats_ref(*(c_end-1)) && --c_end != c_beg);
  if (c_beg == c_end) throw dnmt_error("cigar eats no ref");

  for (auto c_itr = c_beg; c_itr != c_end; ++c_itr)
    if (bam_cigar_op(*c_itr) == BAM_CSOFT_CLIP)
      *c_itr = to_insertion(*c_itr);
}

static void
fix_external_insertion(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 2) return;

  auto c_itr = cigar;
  const auto c_end = c_itr + n_cigar;

  for (; !eats_ref(*c_itr) && c_itr != c_end; ++c_itr)
    *c_itr = to_softclip(*c_itr);

  if (c_itr == c_end) throw dnmt_error("cigar eats no ref");

  c_itr = cigar + n_cigar - 1;
  for (; !eats_ref(*c_itr) && c_itr != cigar; --c_itr)
    *c_itr = to_softclip(*c_itr);
}

static size_t
merge_cigar_ops(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 2) return n_cigar;
  auto c_itr1 = cigar;
  auto c_end = c_itr1 + n_cigar;
  auto c_itr2 = c_itr1 + 1;
  auto op1 = bam_cigar_op(*c_itr1);
  while (c_itr2 != c_end) {
    auto op2 = bam_cigar_op(*c_itr2);
    if (op1 == op2) {
      *c_itr1 = bam_cigar_gen(bam_cigar_oplen(*c_itr1) +
                                  bam_cigar_oplen(*c_itr2), op1);
    }
    else {
      *(++c_itr1) = *c_itr2;
      op1 = op2;
    }
    ++c_itr2;
  }
  // another increment to move past final "active" element for c_itr1
  ++c_itr1;
  return std::distance(cigar, c_itr1);
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
const uint8_t byte_revcom_table[] = {
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



// places seq of b at the end of seq of c
// assumes 0 < c_seq_len - b_seq_len <= a_seq_len
// also assumes that c_seq_len has been figured out
// Also assumes the number of bytes allocated to sequence potion of c->data
// has been set to ceil((a_used_len + b_seq_len) / 2.0) where
// a_used_len = c_seq_len - b_seq_len
static void
merge_by_byte(const bam_rec &a, const bam_rec &b, bam_rec &c) {
  // ADS: (todo) need some functions for int_ceil and is_odd
  const size_t b_seq_len = b.get_l_qseq();
  const size_t c_seq_len = a.get_l_qseq();
  const size_t a_used_len = c_seq_len - b_seq_len;

  const bool is_a_odd = a_used_len % 2 == 1;
  const bool is_b_odd = b_seq_len % 2 == 1;
  const bool is_c_odd = c_seq_len % 2 == 1;

  const size_t a_num_bytes = (a_used_len + 1) / 2;
  const size_t b_num_bytes = (b_seq_len + 1) / 2;

  const size_t b_offset = is_a_odd && is_b_odd;

  const auto a_seq = a.get_qseq_const();
  const auto b_seq = b.get_qseq_const();
  auto c_seq = c.get_qseq();

  std::copy_n(a_seq, a_num_bytes, c_seq);
  if (is_a_odd) {
    // c_seq looks like either [ aa aa aa aa ]
    //                      or [ aa aa aa a- ]
    c_seq[a_num_bytes - 1] &= 0xf0;
    c_seq[a_num_bytes - 1] |= is_b_odd ?
        byte_revcom_table[b_seq[b_num_bytes - 1]] :
        byte_revcom_table[b_seq[b_num_bytes - 1]] >> 4;
  }
  if (is_c_odd) {
    // c_seq looks like either [ aa aa aa aa ]
    //                      or [ aa aa aa ab ]
    for (size_t i = 0; i < b_num_bytes - 1; i++) {
      c_seq[a_num_bytes + i] =
          (byte_revcom_table[b_seq[b_num_bytes - i - 1]] << 4) |
          (byte_revcom_table[b_seq[b_num_bytes - i - 2]] >> 4);
    }
    c_seq[a_num_bytes + b_num_bytes - 1] = byte_revcom_table[b_seq[0]] << 4;
    // Here, c_seq is either [ aa aa aa aa bb bb bb b- ] (a even; b odd)
    //                    or [ aa aa aa ab bb bb bb b- ] (a odd; b odd)
  }
  else {
    for (size_t i = 0; i < b_num_bytes - b_offset; i++) {
      c_seq[a_num_bytes + i] =
          byte_revcom_table[b_seq[b_num_bytes - i - 1 - b_offset]];
    }
    // Here, c_seq is either [ aa aa aa aa bb bb bb bb ] (a even and b even)
    //                    or [ aa aa aa ab bb bb bb    ] (a odd and b odd)
  }
}


bool
check_format_in_header(const string &input_format, const string &inputfile) {
  bam_infile hts(inputfile);
  if (!hts) throw dnmt_error("error opening file: " + inputfile);

  //sam_hdr_t *hdr = sam_hdr_read(hts);
  //if (!hdr) throw dnmt_error("failed to read header: " + inputfile);

  auto begin_hdr = sam_hdr_str(hts.get_hdr_nonconst());
  auto end_hdr = begin_hdr + std::strlen(begin_hdr);
  auto it = std::search(begin_hdr, end_hdr,
                        begin(input_format), end(input_format),
                        [](const unsigned char a, const unsigned char b) {
                          return std::toupper(a) == std::toupper(b);
                        });

  return it != end_hdr;
}





int
merge_non_overlap(const bam_rec &a, const bam_rec &b,
                  const uint32_t spacer, bam_rec &c) {
  /* make the cigar string */
  // collect info about the cigar strings
  const uint32_t *a_cig = a.get_cigar();
  const uint32_t a_ops = a.get_n_cigar();
  const uint32_t *b_cig = b.get_cigar();
  const uint32_t b_ops = b.get_n_cigar();
  // allocate the new cigar string
  const uint32_t c_ops = a_ops + b_ops + 1;
  uint32_t *c_cig = (uint32_t *)calloc(c_ops, sizeof(uint32_t));
  // concatenate the new cigar strings with a "skip" in the middle
  memcpy(c_cig, a_cig, a_ops * sizeof(uint32_t));
  c_cig[a_ops] = bam_cigar_gen(spacer, BAM_CREF_SKIP);
  memcpy(c_cig + a_ops + 1, b_cig, b_ops * sizeof(uint32_t));
  /* done with cigars */

  const size_t a_seq_len = a.get_l_qseq();
  const size_t b_seq_len = b.get_l_qseq();
  const size_t c_seq_len = a_seq_len + b_seq_len;

  // get the template length from the cigar
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag: only need to keep strand and single-end info
  const uint16_t flag = a.get_flag() & (BAM_FREAD1 |
                                        BAM_FREAD2 |
                                        BAM_FREVERSE);

  int ret =
      bam_set1_wrapper(c.get_ptr(),
                       a.get_m_qname(),
                       a.qname_ptr(),
                       flag, // flags (no PE; revcomp info)
                       a.get_tid(),
                       a.get_pos(),
                       a.get_qual(), // mapq from a (consider update)
                       c_ops,        // merged cigar ops
                       c_cig,        // merged cigar
                       -1,           // (no mate)
                       -1,           // (no mate)
                       isize,        // TLEN (relative to reference; SAM docs)
                       c_seq_len,    // merged sequence length
                       8);           // enough for 2 tags of 1 byte value?
  free(c_cig);
  if (ret < 0) throw dnmt_error(ret, "bam_set1 in merge_non_overlap");

  merge_by_byte(a, b, c);

  /* add the tags */
  const int64_t nm = (bam_aux2i(bam_aux_get(a.get_ptr_const(), "NM")) +
                      bam_aux2i(bam_aux_get(b.get_ptr_const(), "NM")));
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c.get_ptr(), "NM", nm);
  if (ret < 0) throw dnmt_error(ret, "merge_non_overlap:bam_aux_update_int");

  const uint8_t cv = bam_aux2A(bam_aux_get(a.get_ptr_const(), "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c.get_ptr(), "CV", 'A', 1, &cv);
  if (ret < 0) throw dnmt_error(ret, "merge_non_overlap:bam_aux_append");

  return ret;
}


int
merge_overlap(const bam_rec &a, const bam_rec &b,
              const uint32_t head, bam_rec &c) {
  assert(head > 0);

  const uint32_t *a_cig = a.get_cigar();
  const uint32_t a_ops = a.get_n_cigar();

  const uint32_t *b_cig = b.get_cigar();
  const uint32_t b_ops = b.get_n_cigar();

  uint32_t part_op = 0;
  uint32_t c_cur = get_full_and_partial_ops(a_cig, a_ops, head, &part_op);
  // ADS: hack here because the get_full_and_partial_ops doesn't do
  // exactly what is needed for this.
  const bool use_partial = (c_cur < a_ops && part_op > 0);

  // check if the middle op would be the same
  const bool merge_mid =
      (use_partial > 0 ?
      bam_cigar_op(a_cig[c_cur]) == bam_cigar_op(b_cig[0]) :
      bam_cigar_op(a_cig[c_cur - 1]) == bam_cigar_op(b_cig[0]));

  // c_ops: include the prefix of a_cig we need; then add for the
  // partial op; subtract for the identical op in the middle; finally
  // add the rest of b_cig.
  const uint32_t c_ops = c_cur + use_partial - merge_mid + b_ops;

  uint32_t *c_cig = (uint32_t *)calloc(c_ops, sizeof(uint32_t));
  // std::fill(c_cig, c_cig + c_ops, std::numeric_limits<uint32_t>::max());
  memcpy(c_cig, a_cig, c_cur * sizeof(uint32_t));

  if (use_partial) {
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
    c_cur++; // index of dest for copying b_cig; faciltates corner case
  }
  // Here we get the length of a's sequence part contribution to c's
  // sequence before the possibility of merging the last entry with
  // the first entry in b's cigar. This is done with the cigar, so
  // everything depends on the "use_partial"
  const size_t a_seq_len = bam_cigar2qlen(c_cur, c_cig);
  /* ADS: above the return type of bam_cigar2qlen is uint64_t, but
     according to the source as of 05/2023 it cannot become
     negative; no possible error code returned */

  if (merge_mid) // update the middle op if it's the same
    c_cig[c_cur - 1] = bam_cigar_gen(bam_cigar_oplen(c_cig[c_cur - 1]) +
                                     bam_cigar_oplen(b_cig[0]),
                                     bam_cigar_op(b_cig[0]));
  // copy the cigar from b into c
  memcpy(c_cig + c_cur, b_cig + merge_mid,
         (b_ops - merge_mid)*sizeof(uint32_t));
  /* done with cigar string here */

  const uint32_t c_seq_len = a_seq_len + b.get_l_qseq();

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag only needs to worry about strand and single-end stuff
  uint16_t flag = (a.get_flag() & (BAM_FREAD1 |
                                   BAM_FREAD2 |
                                   BAM_FREVERSE));

  int ret = bam_set1_wrapper(c.get_ptr(),
                             a.get_m_qname(),
                             a.qname_ptr(),
                             flag, // (no PE; revcomp info)
                             a.get_tid(),
                             a.get_pos(),
                             a.get_qual(), // mapq from "a" (consider update)
                             c_ops,        // merged cigar ops
                             c_cig,        // merged cigar
                             -1,           // (no mate)
                             -1,           // (no mate)
                             isize,        // updated
                             c_seq_len,    // merged sequence length
                             8);           // enough for 2 tags?
  free(c_cig);
  if (ret < 0) throw dnmt_error(ret, "bam_set1_wrapper in merge_overlap");
  // Merge the sequences by bytes
  merge_by_byte(a, b, c);

  // add the tag for mismatches
  const int64_t nm = (bam_aux2i(bam_aux_get(a.get_ptr_const(), "NM")) +
                      bam_aux2i(bam_aux_get(b.get_ptr_const(), "NM")));
  ret = bam_aux_update_int(c.get_ptr(), "NM", nm);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_update_int in merge_overlap");

  // add the tag for conversion
  const uint8_t cv = bam_aux2A(bam_aux_get(a.get_ptr_const(), "CV"));
  ret = bam_aux_append(c.get_ptr(), "CV", 'A', 1, &cv);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_append in merge_overlap");

  return ret;
}




int
truncate_overlap(const bam_rec &a, const uint32_t overlap, bam_rec &c) {

  const uint32_t *a_cig = a.get_cigar();
  const uint32_t a_ops = a.get_n_cigar();

  uint32_t part_op = 0;
  const uint32_t c_cur =
    get_full_and_partial_ops(a_cig, a_ops, overlap, &part_op);

  // ADS: hack here because the get_full_and_partial_ops doesn't do
  // exactly what is needed for this.
  const bool use_partial = (c_cur < a_ops && part_op > 0);

  const uint32_t c_ops = c_cur + use_partial;
  uint32_t *c_cig = (uint32_t *)calloc(c_ops, sizeof(uint32_t));

  // ADS: replace this with a std::copy
  memcpy(c_cig, a_cig, c_cur * sizeof(uint32_t));
  // ADS: warning, if !use_partial, the amount of part_op used below
  // would make no sense.
  if (use_partial)
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
  /* after this point the cigar is set and should decide everything */

  const uint32_t c_seq_len = bam_cigar2qlen(c_ops, c_cig);
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag only needs to worry about strand and single-end stuff
  const uint16_t flag =
      a.get_flag() & (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);

  int ret = bam_set1_wrapper(c.get_ptr(),
                             a.get_m_qname(),
                             a.qname_ptr(),
                             flag, // flags (SR and revcomp info)
                             a.get_tid(),
                             a.get_pos(),
                             a.get_qual(),
                             c_ops,     // merged cigar ops
                             c_cig,     // merged cigar
                             -1,        // (no mate)
                             -1,        // (no mate)
                             isize,     // rlen from new cigar
                             c_seq_len, // truncated seq length
                             8);        // enough for the 2 tags?
  if (ret < 0) throw dnmt_error(ret, "bam_set1_wrapper");
  // ADS: might it be better to fill `c->data` directly?
  free(c_cig);

  auto c_seq = c.get_qseq();
  size_t num_bytes_to_copy = (c_seq_len + 1)/2;
  std::copy_n(bam_get_seq(a.get_ptr_const()), num_bytes_to_copy, c_seq);

  /* add the tags */
  const int64_t nm = bam_aux2i(bam_aux_get(a.get_ptr_const(), "NM")); // ADS: do better here!
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c.get_ptr(), "NM", nm);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_update_int");

  const uint8_t conversion = bam_aux2A(bam_aux_get(a.get_ptr_const(), "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c.get_ptr(), "CV", 'A', 1, &conversion);
  if (ret < 0) throw dnmt_error(ret, "bam_aux_append");

  return ret;
}

size_t
correct_cigar(bam_rec &b) {
  /* This function will change external insertions into soft clip
     operations. Not sure why those would be present. It will also
     change internal soft-clip operations into insertions. This could
     be needed if soft-clipped ends of reads were moved to the middle
     of a merged fragment. Finally, it will collapse adjacent
     identical operations. None of this impacts the seq/qual/aux which
     get moved as a block */

  uint32_t *cigar = b.get_cigar();
  size_t n_cigar = b.get_n_cigar();
  fix_external_insertion(n_cigar, cigar);
  fix_internal_softclip(n_cigar, cigar);

  // merge identical adjacent cigar ops and get new number of ops
  n_cigar = merge_cigar_ops(n_cigar, cigar);
  // difference in bytes to shift the internal data
  const size_t delta = (b.get_n_cigar() - n_cigar) * sizeof(uint32_t);
  if (delta > 0) { // if there is a difference; do the shift
    auto data_end = b.get_data_end();
    std::copy(b.get_qseq_const(), data_end, b.get_qseq() - delta);
    b.get_ptr()->core.n_cigar = n_cigar; // and update number of cigar ops
  }
  return delta;
}

int
keep_better_end(const bam_rec &a, const bam_rec &b, bam_rec &c) {
  if (a.rlen_from_cigar() >= b.rlen_from_cigar()) 
    c.get_ptr_ref() = bam_copy1(c.get_ptr(), a.get_ptr_const());
  else
    c.get_ptr_ref() = bam_copy1(c.get_ptr(), b.get_ptr_const());
  c.get_ptr()->core.mtid = -1;
  c.get_ptr()->core.mpos = -1;
  c.get_ptr()->core.isize = c.rlen_from_cigar();
  c.get_ptr()->core.flag &= (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);
  return 0;
}

