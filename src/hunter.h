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

#ifndef HUNTER_H
#define HUNTER_H

#include <sdsl/suffix_arrays.hpp>

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
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

#include <nlohmann/json.hpp>

#include "neighbors.h"
#include "util.h"
#include "needle.h"

using namespace sdsl;

namespace dicey
{

  struct HunterConfig {
    bool indel;
    bool reverse;
    uint32_t distance;
    std::size_t max_locations;
    std::string sequence;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
  };


  struct DnaHit {
    int32_t score;
    uint32_t chr;
    uint32_t start;
    uint32_t end;
    char strand;
    std::string refalign;
    std::string queryalign;
    
    DnaHit(int32_t sc, uint32_t const refIndex, uint32_t const s, uint32_t const e, char const orient, std::string const& ra, std::string const& qa) : score(sc), chr(refIndex), start(s), end(e), strand(orient), refalign(ra), queryalign(qa) {}
  };

  template<typename TRecord>
  struct SortDnaHit : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& a, TRecord const& b) const {
      return ((a.score > b.score) || ((a.score == b.score) && (a.chr < b.chr)) || ((a.score == b.score) && (a.chr == b.chr) && (a.start < b.start)));
    }
  };


  template<typename TScore>
  inline int32_t
  needleScore(std::string const& a, std::string const& b, TScore const& sc) {
    int32_t score = 0;
    for(uint32_t i = 0; ((i<a.size()) && (i<b.size())); ++i) {
      if (a[i] == b[i]) score += sc.match;
      else score += sc.mismatch;
    }
    return score;
  }

  template<typename TConfig>
  inline void
  jsonDnaHitOut(TConfig const& c, std::vector<std::string> const& qn, std::vector<DnaHit> const& ht) {
    // Output file
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.outfile.c_str(), std::ios_base::out | std::ios_base::binary));
    rcfile << '[';
    
    for(uint32_t i = 0; i < ht.size(); ++i) {
      if (i>0) rcfile << ',';
      nlohmann::json j;
      j["score"] = ht[i].score;
      j["chr"] = qn[ht[i].chr];
      j["start"] = ht[i].start;
      j["end"] = ht[i].end;
      j["strand"] = std::string(1, ht[i].strand);
      j["refalign"] = ht[i].refalign;
      j["queryalign"] = ht[i].queryalign;
      rcfile << j.dump();
    }
    rcfile << ']' << std::endl;
    rcfile.pop();
  }
  
  int hunter(int argc, char** argv) {
    HunterConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("hits.json.gz"), "output file")
      ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(1000), "max. number of matches")
      ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
      ("hamming,n", "use hamming neighborhood instead of edit distance")
      ("forward,f", "only forward matches")
       ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<std::string>(&c.sequence), "input DNA sequence")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] -g Danio_rerio.fa.gz CATTACTAACGTTCAGT" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Cmd switches
    if (!vm.count("hamming")) c.indel = true;
    else c.indel = false;
    if (!vm.count("forward")) c.reverse = true;
    else c.reverse = false;

    // Upper case
    c.sequence = boost::to_upper_copy(c.sequence);
    c.sequence = replaceNonDna(c.sequence);
    std::string revSequence = c.sequence;
    reverseComplement(revSequence);

    
    // Check distance
    if (c.distance >= c.sequence.size()) c.distance = c.sequence.size() - 1;

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "dicey ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Parse chromosome lengths
    std::vector<uint32_t> seqlen;
    std::vector<std::string> seqname;
    uint32_t nseq = getSeqLenName(c, seqlen, seqname);
    if (!nseq) {
      std::cerr << "Could not retrieve sequence lengths!" << std::endl;
      return 1;
    }

    // Reference index
    csa_wt<> fm_index;  
    boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
    std::string index_file = op.string() + ".fm9";
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FM-Index" << std::endl;
    if (!load_from_file(fm_index, index_file)) {
      std::cerr << "Index cannot be loaded!" << std::endl;
      return 1;
    }

    // Define alphabet
    typedef std::set<char> TAlphabet;
    char tmp[] = {'A', 'C', 'G', 'T'};
    TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));

    // Generate neighbors
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Generate forward sequence neighborhood " << std::flush;
    typedef std::set<std::string> TStringSet;
    typedef std::vector<TStringSet> TFwdRevSearchSets;
    TFwdRevSearchSets fwrv(2, TStringSet());
    neighbors(c.sequence, alphabet, c.distance, c.indel, fwrv[0]);
    std::cout << "(#n=" << fwrv[0].size() << ")" << std::endl;
    // Debug
    for(TStringSet::iterator it = fwrv[0].begin(); it != fwrv[0].end(); ++it) std::cerr << *it << std::endl;

    // Reverse complement set?
    if (c.reverse) {
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Generate reverse sequence neighborhood " << std::flush;
      neighbors(revSequence, alphabet, c.distance, c.indel, fwrv[1]);
      std::cout << "(#n=" << fwrv[1].size() << ")" << std::endl;
    }

    // Search
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Search genomic matches" << std::endl;
    uint32_t hits = 0;
    std::vector<DnaHit> ht;
    for(uint32_t fwrvidx = 0; fwrvidx < fwrv.size(); ++fwrvidx) {
      for(typename TStringSet::const_iterator it = fwrv[fwrvidx].begin(); ((it != fwrv[fwrvidx].end()) && (hits < c.max_locations)); ++it) {
	std::string query = *it;
	std::size_t m = query.size();
	std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
	if (occs > 0) {
	  auto locations = locate(fm_index, query.begin(), query.begin() + m);
	  std::sort(locations.begin(), locations.end());
	  for(std::size_t i = 0; ((i < std::min(occs, c.max_locations)) && (hits < c.max_locations)); ++i) {
	    int64_t bestPos = locations[i];
	    int64_t cumsum = 0;
	    uint32_t refIndex = 0;
	    for(; bestPos >= cumsum + seqlen[refIndex]; ++refIndex) cumsum += seqlen[refIndex];
	    uint32_t chrpos = bestPos - cumsum;
	    auto s = extract(fm_index, locations[i], locations[i] + m -1);

	    // Create DNA hit
	    std::string genomicseq = s.substr(0, m);
	    DnaScore<int32_t> sc(1, -1, -1, -1);
	    typedef boost::multi_array<char, 2> TAlign;
	    AlignConfig<false, false> global;
	    if (fwrvidx == 0) {
	      if (c.indel) {
		TAlign align;
		int32_t score = needle(genomicseq, c.sequence , align, global, sc);
		std::string refalign = "";
		std::string queryalign = "";
		for(uint32_t j = 0; j<align.shape()[1]; ++j) {
		  refalign += align[0][j];
		  queryalign += align[1][j];
		}
		ht.push_back(DnaHit(score, refIndex, chrpos+1, chrpos+query.size(), '+', refalign, queryalign));
	      } else {
		int32_t score = needleScore(genomicseq, c.sequence, sc);
		ht.push_back(DnaHit(score, refIndex, chrpos+1, chrpos+query.size(), '+', genomicseq, c.sequence));
	      }
	    } else {
	      if (c.indel) {
		TAlign align;
		int32_t score = needle(genomicseq, revSequence, align, global, sc);
		std::string refalign = "";
		std::string queryalign = "";
		for(uint32_t j = 0; j<align.shape()[1]; ++j) {
		  refalign += align[0][j];
		  queryalign += align[1][j];
		}
		ht.push_back(DnaHit(score, refIndex, chrpos+1, chrpos+query.size(), '-', refalign, queryalign));
	      } else {
		int32_t score = needleScore(genomicseq, revSequence, sc);
		ht.push_back(DnaHit(score, refIndex, chrpos+1, chrpos+query.size(), '-', genomicseq, revSequence));
	      }
	    }
	    ++hits;
	  }
	}
      }
    }

    // Sort
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Sort matches" << std::endl;
    std::sort(ht.begin(), ht.end(), SortDnaHit<DnaHit>());

    // Output
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output matches" << std::endl;
    jsonDnaHitOut(c, seqname, ht);
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    
    return 0;
  }

}

#endif
