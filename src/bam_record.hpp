#ifndef BAM_RECORD_HPP
#define BAM_RECORD_HPP

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <string>

#include <cassert>
#include <iostream>
#include <vector>
#include <stdexcept>




class bam_header {
public:
  bam_header() { header = sam_hdr_init(); }

  bam_header(const bam_header &hdr) {
    header = sam_hdr_dup(hdr.header);
    error_code = (header == nullptr) ? -1 : 0;
  }


  ~bam_header() {
    if (header != nullptr) sam_hdr_destroy(header);
  }

  // not sure if this will remain, but I don't think it can hurt
  explicit bam_header(sam_hdr_t *const hdr) {
    header = sam_hdr_dup(hdr);
    error_code = (header == nullptr) ? -1 : 0;
  }

  std::string
  tostring() const {
    return std::string(sam_hdr_str(header));
  }

  const int32_t
  n_targets() const {
    return header->n_targets;
  }


  inline std::string
  target_name(const int32_t tid) {
    return header->target_name[tid];
  }

  // MN: This is under construction. Do not use this. 
  int
  add_line(const std::string &type, const std::vector<std::string> &args);

  sam_hdr_t *header;
  // ADS: Not sure we need this as an instance variable. Isn't an
  // error something we always manage locally? Do we ever want a
  // header object to persist in an erroneous state?
  int error_code;
};

static inline std::string
to_string(const bam_header &hdr) {
  return hdr.tostring();
}

class bam_rec {
public:
  bam_rec() { record = bam_init1(); }

  bam_rec(const bam_rec &aln) {
    record = bam_init1();
    record = bam_copy1(record, aln.record);
  }

  bam_rec &
  operator=(bam_rec rhs) {
    std::swap(record, rhs.record);
    return *this;
  }

  explicit bam_rec(const bam1_t *aln) {
    record = bam_init1();
    record = bam_copy1(record, aln);
  }

  static void
  steal(bam1_t *&aln, bam_rec &br) {
    bam_destroy1(br.record); // ADS: not a great way to do this!
    br.record = aln;
    aln = nullptr;
  }

  void
  copy(const bam_rec &aln) {
    record = bam_copy1(record, aln.record);
  }
  

  ~bam_rec() {
    // ADS: is it possible for a bam_rec to have `record == NULL`? Can
    // that state be reached?
    if (record != nullptr) {
      bam_destroy1(record);
      record = nullptr; // ADS: in theory this should not be needed
    }
  }



  size_t 
  qlen() const {
    return record->core.l_qseq;
  }


  hts_pos_t &
  pos() {
    return record->core.pos;
  }

  const hts_pos_t &
  pos() const {
    return record->core.pos;
  }

  int32_t &
  tid() {
    return record->core.tid;
  }

  const int32_t &
  tid() const {
    return record->core.tid;
  }

  inline const uint32_t *
  get_cigar() const {
    return reinterpret_cast<const uint32_t *>(record->data + 
        record->core.l_qname);
  }

  inline uint32_t *
  get_cigar() {
    return reinterpret_cast<uint32_t *>(record->data + 
        record->core.l_qname);
  }

  inline const uint8_t *
  get_seq() const {
    return record->data + (record->core.n_cigar << 2) + record->core.l_qname;
  }

  inline uint8_t *
  get_seq() {
    return record->data + (record->core.n_cigar << 2) + record->core.l_qname;
  }

  inline const uint8_t *
  get_qual() const {
    return (record->data + (record->core.n_cigar << 2) + 
        record->core.l_qname + ((record->core.l_qseq + 1) >> 1));
  }

  const uint8_t *
  get_aux() const {
    return get_qual() + record->core.l_qseq;
  }

  const size_t
  l_aux() const {
    return bam_get_l_aux(record);
  } 

  inline std::string
  qname() const {
    return bam_get_qname(record);
  }

  size_t
  qlen_from_cigar() const {
    return bam_cigar2qlen(record->core.n_cigar, get_cigar());
  }

  size_t
  rlen_from_cigar() const {
    return bam_cigar2rlen(record->core.n_cigar, get_cigar());
  }

  hts_pos_t
  endpos() const {
    return bam_endpos(record);
  }

  int
  aux_update_int(const std::string &tag, const int64_t val) {
    return bam_aux_update_int(record, tag.substr(0,2).c_str(), val);
  }

  bool
  is_rev() const;

  void
  write_cigar(std::stringstream &ss) const;

  void
  write_mtid(std::stringstream &ss) const;

  void
  write_seq(std::stringstream &ss) const;

  void
  write_qual(std::stringstream &ss) const;

  void
  write_aux(std::stringstream &ss) const;

  std::string
  tostring() const;

  bam1_t *record;
};


inline void
swap(bam_rec &a, bam_rec &b) {
  std::swap(a.record, b.record);
}

template<class T> T &
operator<<(T &out, const bam_rec &br) {
  out << br.tostring(); // there can be an issue returning this
                        // directly; for example if the `class T`
                        // returns a bool from its operator::<<
  return out;
}

class bam_infile {
public:
  bam_infile(): file(nullptr), error_code(0) {}

  ~bam_infile() {
    // ADS: condition below ensures the hdr will be destroyed when the
    // file is closed.
    assert(file->bam_header->ref_count == 0);
    if (file != nullptr) hts_close(file);
  }

  bam_infile(const std::string &filename) {
    file = hts_open(filename.c_str(), "r");
    error_code = (file) ? 0 : -1;
    if (!error_code) {
      fmt = hts_get_format(file);
      file->bam_header = bam_hdr_read(file->fp.bgzf);
      if (!file->bam_header) {
        if (file != nullptr) hts_close(file);
        throw std::runtime_error("Null or invalid header.");
      }
    }
    else {
      if (file != nullptr) hts_close(file);
      throw std::runtime_error("File could not be opened.");
    }
  }

  bam_infile &
  get_bam_rec(bam_rec &br) {
    // ADS: here we can probably assume this bam_file's `file` has
    // htsFile::format == htsFormatCategory::sequence_data
    if (error_code) return *this;

    bam1_t *aln = bam_init1();
    error_code = sam_read1(file, file->bam_header, aln);
    if (error_code < 0)
      bam_destroy1(aln);
    else {
      error_code = 0;
      bam_rec::steal(aln, br);
    }
    return *this;
  }

  operator bool() const { return (error_code == 0); }

  bool
  is_bam_or_sam();

  htsFile *file; // ADS: probably needs a better name
  const htsFormat *fmt;
  int error_code; // ADS: need to define this better
};

class bam_outfile {
public:
  bam_outfile(): file(nullptr), error_code(0), sam(true) {}

  ~bam_outfile() {
    // ADS: condition below ensures the hdr will be destroyed when the
    // file is closed.
    // assert(file->bam_header->ref_count == 0);
    if (file != nullptr) hts_close(file);
  }

  bam_outfile(const std::string &filename, const bam_header &bh, 
      const bool _sam = true) : sam(_sam) {

    if (sam) file = hts_open(filename.c_str(), "w");
    else file = hts_open(filename.c_str(), "bw");
    error_code = (file) ? 0 : -1;
    if (!bh.header) error_code = -1;
    if (!error_code) {
      file->bam_header = sam_hdr_dup(bh.header);
      int tmp = sam_hdr_write(file, file->bam_header);
      if (tmp < 0) error_code = tmp;
    }
  }

  bam_outfile &
  put_bam_rec(const bam_rec &br) {
    // ADS: here we can probably assume this bam_file's `file` has
    // htsFile::format == htsFormatCategory::sequence_data
    if (error_code) {
      std::cout << "Error: " << error_code << std::endl;
      return *this;
    }
    assert(file->bam_header != nullptr);
    int tmp = sam_write1(file, file->bam_header, br.record);
    if (tmp < 0) error_code = tmp;
    return *this;
  }

  operator bool() const { return (error_code == 0); }

  htsFile *file;  // ADS: probably needs a better name
  int error_code; // ADS: need to define this better
  bool sam;
};

template<class T> inline T &
operator<<(T &out, const bam_header &hdr) {
  out << hdr.tostring(); // there can be an issue returning this
                         // directly; for example if the `class T`
                         // returns a bool from its operator::<<
  return out;
}

inline bam_infile &
operator>>(bam_infile &in, bam_rec &br) {
  return in.get_bam_rec(br);
}

inline bam_outfile &
operator<<(bam_outfile &out, const bam_rec &br) {
  return out.put_bam_rec(br);
}


class bam_tpool{
public:
  bam_tpool(const int n, const int qsize) {
    tpool.pool = hts_tpool_init(n);
    tpool.qsize = qsize;
  }

  ~bam_tpool() {
    hts_tpool_destroy(tpool.pool);
  }

  template<class T> int
  set_io(const T &bam_io_file) {
    return(hts_set_thread_pool(bam_io_file.file, &tpool));
  }

  htsThreadPool tpool;
};




#endif
