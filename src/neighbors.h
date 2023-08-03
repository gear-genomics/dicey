#ifndef NEIGHBORS_H
#define NEIGHBORS_H

#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <sdsl/suffix_arrays.hpp>
#include <htslib/faidx.h>

using namespace sdsl;

namespace dicey {

// Special insert function that makes sure strset has no superstrings
template<typename TStrSet>
inline void
_insert(TStrSet& strset, std::string const& s) {
  bool insertS = true;
  for (typename TStrSet::iterator it = strset.begin(); it != strset.end(); ) {
    if (it->find(s) != std::string::npos) strset.erase(it++); // s is a substring of *it, erase *it
    else {
      if (s.find(*it) != std::string::npos) insertS = false; // *it is a substring of s, do not insert s
      ++it;
    }
  }
  if (insertS) strset.insert(s);
}

template<typename TAlphabet, typename TStringSet>
inline void
_neighbors(std::string const& query, TAlphabet const& alphabet, int32_t const inputdist, int32_t dist, bool indel, int32_t pos, uint32_t maxsize, TStringSet& strset) {
  if (strset.size() < maxsize) {
    if (pos < (int32_t) query.size()) {
      if (dist > 0) {
	if (indel) {
	  // Deletion
	  std::string newst = query.substr(0, pos) + query.substr(pos + 1);
	  _neighbors(newst, alphabet, inputdist, dist - 1, indel, pos, maxsize, strset);
	}
      }
      
      // No change, move to next pos
      _neighbors(query, alphabet, inputdist, dist, indel, pos+1, maxsize, strset);
      
      if (dist > 0) {
	// Switch nucleotide
	for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
	  if (*ait != query[pos]) {
	    std::string newst(query);
	    newst[pos] = *ait;
	    _neighbors(newst, alphabet, inputdist, dist - 1, indel, pos+1, maxsize, strset);
	  }
	}
	
	if (indel) {    
	  // Insertion
	  for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
	    std::string ins("N");
	    ins[0] = *ait;
	    std::string newst = query.substr(0, pos) + ins + query.substr(pos);
	    _neighbors(newst, alphabet, inputdist, dist - 1, indel, pos + 1, maxsize, strset);
	  }
	}
      }
    } else {
      // Only insert true neighbors
      if (dist < inputdist) _insert(strset, query);
    }
  }
  if (inputdist == 0) _insert(strset, query);
}
      

template<typename TAlphabet, typename TStringSet>
inline void
neighbors(std::string const& query, TAlphabet const& alphabet, int32_t dist, bool indel, uint32_t maxsize, TStringSet& strset) {
  _neighbors(query, alphabet, dist, dist, indel, 0, maxsize, strset);
}


}

#endif
