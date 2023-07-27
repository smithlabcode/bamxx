/* MIT License
 *
 * Copyright (c) 2023 Andrew Smith and Masaru Nakajima
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef BAM_RECORD_HPP
#define BAM_RECORD_HPP

#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/sam.h>  // from HTSlib

#include "bam_header.hpp"

class bam_rec {
public:
  bam_rec() { record = bam_init1(); }

  bam_rec(const bam_rec &aln) {
    record = bam_init1();
    record = bam_copy1(record, aln.record);
  }

  bam_rec &operator=(bam_rec rhs) {
    std::swap(record, rhs.record);
    return *this;
  }

  explicit bam_rec(const bam1_t *aln) {
    record = bam_init1();
    record = bam_copy1(record, aln);
  }

  ~bam_rec() {
    if (record != nullptr) {
      bam_destroy1(record);
      record = nullptr;
    }
  }

  // Function to check that the core fields are consistent with each
  // other and the data array.
  bool validate(const sam_hdr_t *const hdr) const;

  // length of the query sequence
  const int32_t get_l_qseq() const { return record->core.l_qseq; }

  const int32_t get_qname_length() const {
    // number of letters in the query's name
    // "+ 1" below: extranul counts *extras*; we don't want *any* nulls
    return record->core.l_qname - (record->core.l_extranul + 1);
  }

  void set_seq(const std::string &s);
  void set_seq(const std::string &s, const std::string &q);
  void set_cigar(const size_t updated_n_cigar, const uint32_t *c);

  const bam1_t &get_rec() const { return *record; }

  //MN: temporary accessor
  bam1_t *get_ptr() { return record; }
  //MN: temporary accessor
  bam1_t *&get_ptr_ref() { return record; }
  //MN: temporary accessor
  const bam1_t *get_ptr_const() const { return record; }
  //MN: temporary accessor
  const uint8_t *get_data_end() const {
    return bam_get_aux(record) + bam_get_l_aux(record);
  }

  int l_data() const { return record->l_data; }

  hts_pos_t get_pos() const { return record->core.pos; }

  hts_pos_t get_mpos() const { return record->core.mpos; }

  uint32_t get_n_cigar() const { return record->core.n_cigar; }

  uint16_t get_flag() const {return record->core.flag; }

  uint16_t get_l_qname() const {return record->core.l_qname;}

  uint16_t get_qual() const {return record->core.qual;}

  // Get the pointer to the cigar portion of the record.
  uint32_t *get_cigar() {
    // reinterpret_cast because data is stored as array of uint8_t,
    // but `cigar` part is aligned at 4-bytes and has 4-byte encoding
    return reinterpret_cast<uint32_t *>(record->data + record->core.l_qname);
  }

  const uint32_t *get_cigar() const {
    return reinterpret_cast<uint32_t *>(record->data + record->core.l_qname);
  }


  // ADS: not sure we should have this function, as we don't want to
  // change this without actually changing the cigar itself.

  // Get the pointer to the query sequence part of the record.
  uint8_t *get_qseq() {
    return record->data + (record->core.n_cigar << 2) + record->core.l_qname;
  }

  const uint8_t *get_qseq_const() const {
    return record->data + (record->core.n_cigar << 2) + record->core.l_qname;
  }

  // ADS: likely should be a private function. We likely don't want to
  // allow users to change anything inside there without using a
  // function of this class that operates on the `seq` part and also
  // the parts that come after.

  // get length of auxiliary data
  const size_t l_aux() const { return bam_get_l_aux(record); }

  std::string get_qname() const { return bam_get_qname(record); }

  std::string qname() const { return bam_get_qname(record); }
  char *qname_ptr() const { return bam_get_qname(record); }

  // ADS: use cigar to get `qlen`: length of query sequence
  size_t qlen_from_cigar() const {
    return bam_cigar2qlen(record->core.n_cigar, get_cigar());
  }

  // ADS: use cigar to get `rlen`: range in target covered by query
  size_t rlen_from_cigar() const {
    return bam_cigar2rlen(record->core.n_cigar, get_cigar());
  }

  // ADS: seems like this is `pos + rlen`??
  hts_pos_t get_endpos() const { return bam_endpos(record); }

  /* AUX FIELDS */
  template<typename T> void aux_update_int(const char tag[2], const T val) {
    const int ret = bam_aux_update_int(record, tag, static_cast<int64_t>(val));
    if (ret != 0) throw std::runtime_error("fail aux_update_int");
  }

  template<typename T> void aux_update_float(const char tag[2], const T val) {
    const int ret = bam_aux_update_float(record, tag, static_cast<float>(val));
    if (ret != 0) throw std::runtime_error("fail aux_update_float");
  }

  template<typename T>
  void aux_update_array(const char tag[2], const uint8_t tag_type,
                        const uint32_t items, const T *data) {
    // warning: bam_aux_update_array does not take const data
    const int ret = bam_aux_update_array(record, tag, tag_type, items,
                                         reinterpret_cast<void *>(data));
    if (ret != 0) throw std::runtime_error("fail aux_update_array");
  }

  void aux_update_str(const char tag[2], const uint8_t tag_type,
                      const std::string &data);

  template<typename T_itr> void append_aux(const char tag[2],
                                           const uint8_t aux_type, T_itr b_dat,
                                           T_itr e_dat) {
    // ADS: this doesn't work at all; what should aux_type be?
    const size_t sz = std::distance(b_dat, e_dat) * sizeof(*b_dat);
    std::vector<uint8_t> data(sz, 0);
    std::copy(reinterpret_cast<uint8_t *>(b_dat),
              reinterpret_cast<uint8_t *>(e_dat), begin(data));
    const int ret = bam_aux_append(record, tag, aux_type, sz, data.data());
    if (ret < 0) throw std::runtime_error("fail append_aux");
  }

  void append_aux(const char tag[2], const char aux_type, const std::string &s);

  template<typename T>
  void append_aux(const char tag[2], const char aux_type, T dat) {
    const int ret = bam_aux_append(record, tag, aux_type, sizeof(T),
                                   reinterpret_cast<uint8_t *>(&dat));
    if (ret < 0) throw std::runtime_error("fail append_aux");
  }

  void remove_aux(const char tag[2]);

  int32_t get_mtid() const { return record->core.mtid; }

  int32_t get_tid() const { return record->core.tid; }


  // indirect accessors
  bool is_rev() const;
  bool is_a_rich() const;
  std::string tostring(const sam_hdr_t *const hdr) const;

  // ADS: non-member function to take the data from a struct and be
  // responsible for it
  static void steal(bam1_t *&aln, bam_rec &br) {
    bam_destroy1(br.record); // ADS: not a great way to do this!
    br.record = aln;
    aln = nullptr;
  }

  void revcomp();
  void flip();
private:
  bam1_t *record;

  uint32_t &n_cigar() { return record->core.n_cigar; }


  // Get the pointer to the aux part of the bam record.
  uint8_t *get_aux() { return get_qual() + record->core.l_qseq; }

  // Get the pointer to the quality score string for the record.
  uint8_t *get_qual() {
    return (record->data + (record->core.n_cigar << 2) + record->core.l_qname +
            ((record->core.l_qseq + 1) >> 1));
  }

  // lvalue references
  int32_t &mtid() { return record->core.mtid; }

  int32_t &tid() { return record->core.tid; }

  hts_pos_t &pos() { return record->core.pos; }

  int &l_data() { return record->l_data; }

  int32_t &l_qseq() { return record->core.l_qseq; }

  hts_pos_t &mpos() { return record->core.mpos; }

  hts_pos_t &isize() { return record->core.isize; }

  uint16_t &flag() {return record->core.flag; }
};

#endif
