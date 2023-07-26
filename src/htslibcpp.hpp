
#ifndef HSTLIBCPP
#define HSTLIBCPP

#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "bam_record_ads.hpp"
#include "bam_infile.hpp"

void
standardize_format(const std::string &input_format, bam_rec &aln);

bool
same_name(const bam_rec &a, const bam_rec &b, const size_t suff_len);

bool
are_mates(const bam_rec &one, const bam_rec &two);

bool
check_input_file(const std::string &infile);

std::vector<std::string>
load_read_names(const std::string &inputfile, const size_t n_reads);

bool
check_format_in_header(const std::string &input_format, const std::string &inputfile);

int
merge_non_overlap(const bam_rec &a, const bam_rec &b,
                  const uint32_t spacer, bam_rec &c);

int
merge_overlap(const bam_rec &a, const bam_rec &b,
              const uint32_t head, bam_rec &c);

int
truncate_overlap(const bam_rec &a, const uint32_t overlap, bam_rec &c);

int
keep_better_end(const bam_rec &a, const bam_rec &b, bam_rec &c);

size_t
correct_cigar(bam_rec &b);

#endif
