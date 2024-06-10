#ifndef GTF_H
#define GTF_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace dicey {

  struct GeneInfo {
    bool pcoding;
    std::string id;
    std::string symbol;
    std::string barcode;
    std::string code;
    
    GeneInfo(bool const p, std::string const& idname, std::string const& sym) : pcoding(p), id(idname), symbol(sym), barcode("NNNNNNNNNNNNNNNNNNNN"), code("000000") {}
  };
  
  
  struct IntervalLabel {
    int32_t start;
    int32_t end;
    char strand;
    int32_t lid;

    explicit IntervalLabel(int32_t s) : start(s), end(s+1), strand('*'), lid(-1) {}
    IntervalLabel(int32_t s, int32_t e, char t, int32_t l) : start(s), end(e), strand(t), lid(l) {}
  };

  struct IntervalLabelId {
    int32_t start;
    int32_t end;
    char strand;
    int32_t lid;
    int32_t eid;

    explicit IntervalLabelId(int32_t s) : start(s), end(s+1), strand('*'), lid(-1), eid(-1) {}
    IntervalLabelId(int32_t s, int32_t e, char t, int32_t l, int32_t i) : start(s), end(e), strand(t), lid(l), eid(i) {}
  };
  
  template<typename TRecord>
  struct SortIntervalLabel : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.lid < s2.lid;
    }
  };
  
  template<typename TRecord>
  struct SortIntervalStart : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.start < s2.start;
    }
  };

  inline void
  _insertInterval(std::vector<IntervalLabel>& cr, int32_t s, int32_t e, char strand, int32_t lid, int32_t) {
    // Uniqueness not necessary because we flatten the interval map
    cr.push_back(IntervalLabel(s, e, strand, lid));
  }

  inline void
  _insertInterval(std::vector<IntervalLabelId>& cr, int32_t s, int32_t e, char strand, int32_t lid, int32_t eid) {
    // Check uniqueness
    bool isUnique = true;
    for(uint32_t i = 0; i < cr.size(); ++i) {
      if ((cr[i].start == s) && (cr[i].end == e) && (cr[i].strand == strand) && (cr[i].lid == lid)) {
	isUnique = false;
	break;
      }
    }
    if (isUnique) cr.push_back(IntervalLabelId(s, e, strand, lid, eid));
  }

  
  template<typename TConfig, typename TGenomicRegions, typename TGeneInfo>
  inline int32_t
  parseGTFAll(TConfig const& c, TGenomicRegions& overlappingRegions, TGeneInfo& geneInfo) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "GTF feature parsing" << std::endl;

    // Map IDs to integer
    typedef std::map<std::string, int32_t> TIdMap;
    TIdMap idMap;

    // Keep track of unique exon IDs
    int32_t eid = 0;

    // Parse GTF
    std::ifstream file;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.gtfFile)) {
      file.open(c.gtfFile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else file.open(c.gtfFile.string().c_str(), std::ios_base::in);
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      if ((gline.size()) && (gline[0] == '#')) continue;
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter==tokens.end()) {
	std::cerr << "Empty line in GTF file!" << std::endl;
	return 0;
      }
      std::string chrName=*tokIter++;
      if (c.nchr.find(chrName) == c.nchr.end()) continue;
      int32_t chrid = c.nchr.find(chrName)->second;      
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted GTF file!" << std::endl;
	return 0;
      }
      ++tokIter;
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted GTF file!" << std::endl;
	return 0;
      }
      std::string ft = *tokIter++;
      if (ft == c.feature) {    // Select exons
	if (tokIter != tokens.end()) {
	  int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	  ++tokIter; // score
	  if (tokIter == tokens.end()) {
	    std::cerr << "Corrupted GTF file!" << std::endl;
	    return 0;
	  }
	  char strand = boost::lexical_cast<char>(*tokIter++);
	  ++tokIter; // frame
	  std::string attr = *tokIter;
	  boost::char_separator<char> sepAttr(";");
	  Tokenizer attrTokens(attr, sepAttr);
	  for(Tokenizer::iterator attrIter = attrTokens.begin(); attrIter != attrTokens.end(); ++attrIter) {
	    std::string keyval = *attrIter;
	    boost::trim(keyval);
	    boost::char_separator<char> sepKeyVal(" ");
	    Tokenizer kvTokens(keyval, sepKeyVal);
	    Tokenizer::iterator kvTokensIt = kvTokens.begin();
	    std::string key = *kvTokensIt++;
	    if (key == c.idname) {     // Select gene_id
	      // Protein-coding exon?
	      bool includeExon = false;
	      for(Tokenizer::iterator arIter = attrTokens.begin(); arIter != attrTokens.end(); ++arIter) {
		std::string kvl = *arIter;
		boost::trim(kvl);
		boost::char_separator<char> sKV2(" ");
		Tokenizer kvT2(kvl, sKV2);
		Tokenizer::iterator kvT2It = kvT2.begin();
		std::string procod = *kvT2It++;
		if (procod == "transcript_biotype") {
		  std::string gbio = *kvT2It;
		  if (gbio.size() >= 3) gbio = gbio.substr(1, gbio.size()-2);
		  if (gbio == "protein_coding") includeExon = true;
		}
	      }
	      std::string ensgene = *kvTokensIt;
	      if (ensgene.size() >= 3) ensgene = ensgene.substr(1, ensgene.size()-2); // Trim off the bloody "
	      if ((includeExon) && ((c.computeAll) || (c.geneset.find(ensgene) != c.geneset.end()))) {
		int32_t idval = geneInfo.size();
		typename TIdMap::const_iterator idIter = idMap.find(ensgene);
		if (idIter == idMap.end()) {
		  idMap.insert(std::make_pair(ensgene, idval));
		  // Protein Coding?
		  bool pCode = false;
		  std::string symbol = "n.a.";
		  for(Tokenizer::iterator arIter = attrTokens.begin(); arIter != attrTokens.end(); ++arIter) {
		    std::string kvl = *arIter;
		    boost::trim(kvl);
		    boost::char_separator<char> sKV2(" ");
		    Tokenizer kvT2(kvl, sKV2);
		    Tokenizer::iterator kvT2It = kvT2.begin();
		    std::string procod = *kvT2It++;
		    if (procod == "gene_biotype") {
		      std::string gbio = *kvT2It;
		      if (gbio.size() >= 3) gbio = gbio.substr(1, gbio.size()-2);
		      if (gbio == "protein_coding") pCode = true;
		    }
		    if (procod == "gene_name") {
		      std::string gbio = *kvT2It;
		      if (gbio.size() >= 3) gbio = gbio.substr(1, gbio.size()-2);
		      symbol = gbio;
		    }
		  }
		  geneInfo.push_back(GeneInfo(pCode, ensgene, symbol));
		} else idval = idIter->second;
		// Convert to 0-based and right-open
		if (start == 0) {
		  std::cerr << "GTF is 1-based format!" << std::endl;
		  return 0;
		}
		if (start > end) {
		  std::cerr << "Feature start is greater than feature end!" << std::endl;
		  return 0;
		}
		_insertInterval(overlappingRegions[chrid], start - 1, end - 1, strand, idval, eid++);
	      }
	    }
	  }
	}
      }
    }
    dataIn.pop();
    if (is_gz(c.gtfFile)) dataIn.pop();
    file.close();
    return geneInfo.size();
  }


  template<typename TConfig, typename TGenomicRegions>
  inline int32_t
  parseGTF(TConfig const& c, TGenomicRegions& gRegions, std::vector<GeneInfo>& geneInfo) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Overlapping intervals for each label
    TGenomicRegions overlappingRegions;
    overlappingRegions.resize(gRegions.size(), TChromosomeRegions());
    parseGTFAll(c, overlappingRegions, geneInfo);
    
    // Make intervals non-overlapping for each label
    for(uint32_t refIndex = 0; refIndex < overlappingRegions.size(); ++refIndex) {
      // Sort by ID
      std::sort(overlappingRegions[refIndex].begin(), overlappingRegions[refIndex].end(), SortIntervalLabel<IntervalLabel>());
      int32_t runningId = -1;
      char runningStrand = '*';
      typedef boost::icl::interval_set<uint32_t> TIdIntervals;
      typedef typename TIdIntervals::interval_type TIVal;
      TIdIntervals idIntervals;
      for(uint32_t i = 0; i < overlappingRegions[refIndex].size(); ++i) {
	if (overlappingRegions[refIndex][i].lid != runningId) {
	  for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) {
	    gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
	  }
	  idIntervals.clear();
	  runningId = overlappingRegions[refIndex][i].lid;
	  runningStrand = overlappingRegions[refIndex][i].strand;
	}
	idIntervals.insert(TIVal::right_open(overlappingRegions[refIndex][i].start, overlappingRegions[refIndex][i].end));
      }
      // Process last id
      for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
    }
    
    return geneInfo.size();
  }

  template<typename TConfig, typename TGenomicRegions>
  inline int32_t
  parseFastaSeqs(TConfig const& c, TGenomicRegions& gRegions, std::vector<GeneInfo>& geneInfo) {
    faidx_t* fai = fai_load(c.infile.string().c_str());
    int32_t runningId = 0;
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      std::string chrName = faidx_iseq(fai, refIndex);
      gRegions[refIndex].push_back(IntervalLabel(0, faidx_seq_len(fai, chrName.c_str()) + 1, '+', runningId++));
      geneInfo.push_back(GeneInfo(true, chrName, chrName));
    }
    fai_destroy(fai);
    return geneInfo.size();
  }
  
}

#endif
