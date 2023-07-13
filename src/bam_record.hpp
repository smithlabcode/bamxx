#ifndef BAM_RECORD_HPP 
#define BAM_RECORD_HPP

#include <string>
#include <htslib/sam.h>

#include <iostream>
#include <cassert>
#include <vector>

const uint8_t cigar_shift = 4;
const uint32_t cigar_mask = 0xf;
const std::string cigar_str = "MIDNSHP=XB??????";
const uint8_t sequence_data_enum = 1;
const uint8_t sam_enum = 3;
const uint8_t bam_enum = 4;



const std::vector<std::string> byte2str = {
  "", "A", "C", "NN", "G", "NN", "NN", "NN",
  "T", "NN", "NN", "NN", "NN", "NN", "NN", "N",
  "A", "AA", "AC", "NN", "AG", "NN", "NN", "NN",
  "AT", "NN", "NN", "NN", "NN", "NN", "NN", "AN",
  "C", "CA", "CC", "NN", "CG", "NN", "NN", "NN",
  "CT", "NN", "NN", "NN", "NN", "NN", "NN", "CN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "G", "GA", "GC", "NN", "GG", "NN", "NN", "NN",
  "GT", "NN", "NN", "NN", "NN", "NN", "NN", "GN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "T", "TA", "TC", "NN", "TG", "NN", "NN", "NN",
  "TT", "NN", "NN", "NN", "NN", "NN", "NN", "TN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "NN", "NN", "NN", "NN", "NN", "NN", "NN", "NN",
  "N", "NA", "NC", "NN", "NG", "NN", "NN", "NN",
  "NT", "NN", "NN", "NN", "NN", "NN", "NN", "NN"
};



class bam_header {
public:

  bam_header() {
    header = sam_hdr_init();
  }

  bam_header(const bam_header& hdr) {
    // header = sam_hdr_init();
    header = hdr.header;
    header->ref_count++; // = hdr.header;
  }

  ~bam_header() {
    if (header != NULL)
      sam_hdr_destroy(header);
  }

  // not sure if this will remain, but I don't think it can hurt
  explicit bam_header(const sam_hdr_t *hdr) {
    header = sam_hdr_dup(hdr);
  }

  std::string tostring() const {
    return std::string(sam_hdr_str(header));
  }

  int32_t& n_targets() {return header->n_targets;}
  const int32_t& n_targets() const {return header->n_targets;}

  int32_t& ignore_sam_err() {return header->ignore_sam_err;}
  const int32_t& ignore_sam_err() const {return header->ignore_sam_err;}

  size_t& l_text() {return header->l_text;}
  const size_t& l_text() const {return header->l_text;}

  uint32_t ref_count() {return header->ref_count;}
  const uint32_t ref_count() const {return header->ref_count;}

  sam_hdr_t* get() {return header;}
  const sam_hdr_t* get() const {return header;}

  inline
  std::string target_name(const int32_t tid) {
    return header->target_name[tid];
  }

  sam_hdr_t* header;
};

static inline
std::string
to_string(const bam_header &hdr) {
  return hdr.tostring();
}

class bam_rec {
public:
  bam_rec() {
    record = bam_init1();
  }

  bam_rec(const bam_rec& aln) {
    record = bam_init1();
    record = bam_copy1(record, aln.record);
  }

  explicit bam_rec(const bam1_t *aln) {
    record = bam_init1();
    record = bam_copy1(record, aln);
  }

  static void
  steal(bam1_t* &aln, bam_rec &br) {
    bam_destroy1(br.record); // ADS: not a great way to do this!
    br.record = aln;
    aln = NULL;
  }

  ~bam_rec() {
    // ADS: is it possible for a bam_rec to have `record == NULL`? Can
    // that state be reached?
    if (record != NULL) bam_destroy1(record);
  }

  std::string tostring() const;

  hts_pos_t& pos() {return record->core.pos;}
  const hts_pos_t& pos() const {return record->core.pos;}

  int32_t& tid() {return record->core.tid;}
  const int32_t& tid() const {return record->core.tid;}

  uint16_t& bin() {return record->core.bin;}
  const uint16_t& bin() const {return record->core.bin;}

  uint8_t& qual() {return record->core.qual;}
  const uint8_t& qual() const {return record->core.qual;}

  uint8_t& l_extranul() {return record->core.l_extranul;}
  const uint8_t& l_extranul() const {return record->core.l_extranul;}

  uint16_t& flag() {return record->core.flag;}
  const uint16_t& flag() const {return record->core.flag;}

  uint16_t& l_qname() {return record->core.l_qname;}
  const uint16_t& l_qname() const {return record->core.l_qname;}

  uint32_t& n_cigar() {return record->core.n_cigar;}
  const uint32_t& n_cigar() const {return record->core.n_cigar;}

  int32_t& l_qseq() {return record->core.l_qseq;}
  const int32_t& l_qseq() const {return record->core.l_qseq;}

  int32_t& mtid() {return record->core.mtid;}
  const int32_t& mtid() const {return record->core.mtid;}

  hts_pos_t& mpos() {return record->core.mpos;}
  const hts_pos_t& mpos() const {return record->core.mpos;}

  hts_pos_t& isize() {return record->core.isize;}
  const hts_pos_t& isize() const {return record->core.isize;}

  uint64_t& id() {return record->id;}
  const uint64_t& id() const {return record->id;}

  int& l_data() {return record->l_data;}
  const int& l_data() const {return record->l_data;}

  uint8_t* data() {return record->data;}
  const uint8_t* data() const {return record->data;}

  uint32_t& m_data() {return record->m_data;}
  const uint32_t& m_data() const {return record->m_data;}

  bam1_t* get() {return record;}
  const bam1_t* get() const {return record;}

  inline const uint32_t*
  get_cigar() const {
    return reinterpret_cast<const uint32_t*>(data() + l_qname());
  }

  inline const uint8_t*
  get_seq() const {
    return data() + (n_cigar()<<2) + l_qname();
  }

  inline const uint8_t*
  get_qual() const {
    return (data() + (n_cigar()<<2) + l_qname() + ((l_qseq() + 1)>>1));
  }

  inline const uint8_t*
  get_aux() const {return get_qual() + l_qseq();}

  inline std::string
  qname() const {return bam_get_qname(record);}
  
  inline size_t
  qlen_from_cigar() const {return bam_cigar2qlen(n_cigar(), get_cigar());}

  bam1_t* record;

};



template <class T> T &
operator<<(T &out, const bam_rec &br) {
  out << br.tostring(); // there can be an issue returning this
                        // directly; for example if the `class T`
                        // returns a bool from its operator::<<
  return out;
}





class bam_infile {
public:
  bam_infile() : file(NULL), error_code(0) {}
  ~bam_infile() {
    // ADS: condition below ensures the hdr will be destroyed when the
    // file is closed.
    assert(file->bam_header->ref_count == 0);
    if (file != NULL)
      hts_close(file);
  }

  bam_infile(const std::string &filename) {
    file = hts_open(filename.c_str(), "r");
    error_code = (file) ? 0 : -1;
    if (!error_code) {
      sam_hdr_t *hdr = bam_hdr_read(file->fp.bgzf);
      fmt = hts_get_format(file);
      if (!hdr)
        error_code = -1;
      else {
        assert(hdr->ref_count == 0);
        file->bam_header = hdr;
      }
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

  operator bool() const {
    return (error_code == 0);
  }

  inline bool
  is_bam_or_sam() {
    return fmt->category == sequence_data_enum && 
      (fmt->format == bam_enum || fmt->format == sam_enum);
  }

  htsFile *file;  // ADS: probably needs a better name
  const htsFormat *fmt;
  int error_code; // ADS: need to define this better
};



class bam_outfile {
public:
  bam_outfile() : file(NULL), error_code(0) {}
  ~bam_outfile() {
    // ADS: condition below ensures the hdr will be destroyed when the
    // file is closed.
    // assert(file->bam_header->ref_count == 0);
    if (file != NULL)
      hts_close(file);
  }
  bam_outfile(const std::string &filename,
              bam_header &bh) {
    // ADS: ??? "bh" non-const as its sam_hdr_t* will need to be
    // attached here and might be reference counted?
    file = hts_open(filename.c_str(), "w");
    error_code = (file) ? 0 : -1;
    if (!bh.header)
      error_code = -1;
    if (!error_code) {
      //std::cout << bh.header->ref_count << std::endl;
      file->bam_header = bh.header;
      // using same "sam_hdr_t" so increase the ref_count
      file->bam_header->ref_count++;
      // ADS: do we care about this: add_pg_line??
      int tmp = sam_hdr_write(file, file->bam_header);
      if (tmp < 0)
        error_code = tmp;
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
    assert(file->bam_header != NULL);
    //std::cout << br.tostring() << std::endl;
    int tmp = sam_write1(file, file->bam_header, br.record);
    if (tmp < 0)
      error_code = tmp;
    return *this;
  }

  operator bool() const {
    return (error_code == 0);
  }

  htsFile* file;  // ADS: probably needs a better name
  int error_code; // ADS: need to define this better
};


template <class T> 
inline T &
operator<<(T &out, const bam_header &hdr) {
  out << hdr.tostring(); // there can be an issue returning this
                         // directly; for example if the `class T`
                         // returns a bool from its operator::<<
  return out;
}

inline
bam_infile &
operator>>(bam_infile &in, bam_rec &br) {
  return in.get_bam_rec(br);
}

inline
bam_outfile &
operator<<(bam_outfile &out, const bam_rec &br) {
  return out.put_bam_rec(br);
}

#endif
