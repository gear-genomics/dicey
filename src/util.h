#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

namespace dicey
{

  inline bool is_gz(boost::filesystem::path const& f) {
    std::ifstream bfile(f.string().c_str(), std::ios_base::binary | std::ios::ate);
    bfile.seekg(0, std::ios::beg);
    char byte1;
    bfile.read(&byte1, 1);
    char byte2;
    bfile.read(&byte2, 1);
    bfile.close();
    if ((byte1 == '\x1F') && (byte2 == '\x8B')) return true;
    else return false;
  }

  inline char
  complement(char n) {
    switch(n) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    }
    return 'N';
  }

  inline void
  revcomplement(std::string& nucs) {
    for(std::string::iterator it = nucs.begin(); it != nucs.end(); ++it) *it = complement(*it);
    std::reverse(nucs.begin(), nucs.end());
  }

  inline double
  gccontent(std::string& nucs) {
    uint32_t gc = 0;
    for(std::string::iterator it = nucs.begin(); it != nucs.end(); ++it) {
      if (*it == 'N') return -1;
      else if ((*it == 'C') || (*it == 'G')) ++gc;
    }
    return (double) gc / (double) nucs.size();
  }
  
  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }
  
  inline bool
  nContent(std::string const& s) {
    for(uint32_t i = 0; i < s.size(); ++i) {
      if ((s[i] == 'N') || (s[i] == 'n')) return true;
    }
    return false;
  }

  template<typename TConfig>
  inline bool
  chrNoData(TConfig const& c, uint32_t const refIndex, hts_idx_t const* idx) {
    // Check we have mapped reads on this chromosome
    std::string suffix("cram");
    std::string str(c.bamFile.string());
    if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) return false;
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
    if (mapped) return false;
    else return true;
  }    
  
  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }
  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }


  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }

  inline uint32_t halfAlignmentLength(bam1_t const* rec) {
    return (alignmentLength(rec) / 2);
  }
  
  template<typename TConfig>
  inline int32_t
  getSeqLenName(TConfig const& c, std::vector<uint32_t>& seqlen, std::vector<std::string>& seqname) {
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Input reference file is missing: " << c.genome.string() << std::endl;
      return 0;
    }
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (fai == NULL) {
      if (fai_build(c.genome.string().c_str()) == -1) {
	std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	return 0;
      } else fai = fai_load(c.genome.string().c_str());
    }
    seqlen.resize(faidx_nseq(fai));
    seqname.resize(faidx_nseq(fai), "");
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      std::string sqn(faidx_iseq(fai, refIndex));
      seqlen[refIndex] = faidx_seq_len(fai, sqn.c_str()) + 1;
      seqname[refIndex] = sqn;
    }
    if (fai != NULL) fai_destroy(fai);
    return seqlen.size();
  }

  inline std::string
  replaceNonDna(std::string const& str, std::vector<std::string>& msg) {
    std::string out;
    for(uint32_t i = 0; i<str.size();++i) {
      if ((str[i] == 'A') || (str[i] == 'C') || (str[i] == 'G') || (str[i] == 'T')) out = out.append(str, i, 1);
      else {
	msg.push_back("Warning: Non-DNA character in nucleotide sequence detected and replaced by 'N'!");
	out = out.append("N");
      }
    }
    return out;
  }


}

#endif
