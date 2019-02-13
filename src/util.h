/*
============================================================================
Dicey: In-silico PCR
============================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <htslib/faidx.h>

namespace dicey
{

  inline bool
  nContent(std::string const& s) {
    for(uint32_t i = 0; i < s.size(); ++i) {
      if ((s[i] == 'N') || (s[i] == 'n')) return true;
    }
    return false;
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
