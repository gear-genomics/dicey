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

#include <sdsl/suffix_arrays.hpp>
#include <htslib/faidx.h>

using namespace sdsl;

namespace dicey {

template<typename TStrSet>
inline void
_insert(TStrSet& strset, std::string const& s, bool indel) {
  if (!indel) {
    strset.insert(s);
    return;
  }
  bool insertS = true;
  for (typename TStrSet::iterator it = strset.begin(); it != strset.end(); ) {
    if (it->find(s) != std::string::npos) strset.erase(it++);
    else {
      if (s.find(*it) != std::string::npos) insertS = false;
      ++it;
    }
  }
  if (insertS) strset.insert(s);
}

template<typename TAlphabet, typename TStringSet>
inline void
_neighbors(std::string& query, TAlphabet const& alphabet, int32_t const inputdist, int32_t dist, bool indel, int32_t pos, uint32_t maxsize, TStringSet& strset) {
  if (strset.size() >= maxsize) return;
  if (pos < (int32_t) query.size()) {
    if ((dist > 0) && (indel)) {
      // Deletion
      std::string newst = query.substr(0, pos) + query.substr(pos + 1);
      _neighbors(newst, alphabet, inputdist, dist - 1, indel, pos, maxsize, strset);
    }

    // No change
    _neighbors(query, alphabet, inputdist, dist, indel, pos+1, maxsize, strset);

    if (dist > 0) {
      char orig = query[pos];
      for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
	if (*ait != orig) {
	  query[pos] = *ait;
	  _neighbors(query, alphabet, inputdist, dist - 1, indel, pos+1, maxsize, strset);
	}
      }
      query[pos] = orig;

      if (indel) {
	// Insertion
	for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
	  std::string newst = query.substr(0, pos) + std::string(1, *ait) + query.substr(pos);
	  _neighbors(newst, alphabet, inputdist, dist - 1, indel, pos + 1, maxsize, strset);
	}
      }
    }
  } else {
    // Only insert true neighbors
    if (dist < inputdist) _insert(strset, query, indel);
  }
}


template<typename TAlphabet, typename TStringSet>
inline void
neighbors(std::string const& query, TAlphabet const& alphabet, int32_t dist, bool indel, uint32_t maxsize, TStringSet& strset) {
  std::string q(query);
  _insert(strset, q, indel);
  _neighbors(q, alphabet, dist, dist, indel, 0, maxsize, strset);
}


}

#endif
