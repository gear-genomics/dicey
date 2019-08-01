#ifndef BED_H
#define BED_H

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include <boost/progress.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>


namespace dicey
{


  template<typename TRegionsGenome>
  inline int32_t
  _parseBedIntervals(std::string const& filename, std::vector<std::string> const& chrName, TRegionsGenome& bedRegions) {
    typedef typename TRegionsGenome::value_type TChrIntervals;
    typedef typename TChrIntervals::interval_type TIVal;

    int32_t intervals = 0;
    bedRegions.resize(chrName.size(), TChrIntervals());
    std::ifstream chrFile(filename.c_str(), std::ifstream::in);
    if (chrFile.is_open()) {
      while (chrFile.good()) {
	std::string chrFromFile;
	getline(chrFile, chrFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(chrFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrom = *tokIter++;
	  int32_t tid = -1;
	  for(uint32_t refIndex = 0; refIndex < chrName.size(); ++refIndex) {
	    if (chrom == chrName[refIndex]) {
	      tid = refIndex;
	      break;
	    }
	  }
	  if (tid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      bedRegions[tid].insert(TIVal::right_open(start, end));
	      ++intervals;
	    }
	  }
	}
      }
      chrFile.close();
    }
    return intervals;
  }


}

#endif
