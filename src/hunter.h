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
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
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
    bool hasOutfile;
    uint32_t distance;
    uint32_t maxNeighborhood;
    std::size_t pre_context;
    std::size_t post_context;
    std::size_t max_locations;
    std::string sequence;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
  };


  struct DnaHit {
    int32_t score;
    uint32_t chr;
    uint32_t start;
    char strand;
    std::string refalign;
    std::string queryalign;
    
    DnaHit(int32_t sc, uint32_t const refIndex, uint32_t const s, char const orient, std::string const& ra, std::string const& qa) : score(sc), chr(refIndex), start(s), strand(orient), refalign(ra), queryalign(qa) {}
  };

  template<typename TRecord>
  struct SortDnaHit : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& a, TRecord const& b) const {
      return ((a.score > b.score) || ((a.score == b.score) && (a.chr < b.chr)) || ((a.score == b.score) && (a.chr == b.chr) && (a.start < b.start)));
    }
  };


  template<typename TAlign>
  inline uint32_t
  _trailGap(TAlign const& align) {
    uint32_t lastAlignedPos = align.shape()[1] - 1;
    for(uint32_t j = 0; j<align.shape()[1]; ++j) {
      if (align[1][j] != '-') lastAlignedPos = j;
    }
    return (align.shape()[1] - lastAlignedPos - 1);
  }

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

  inline uint32_t
  _nucleotideLength(std::string const& seq) {
    uint32_t s = 0;
    for(uint32_t i = 0; i < seq.size(); ++i) {
      if (seq[i] != '-') ++s;
    }
    return s;
  }
  
  template<typename TConfig, typename TStream>
  inline void
  writeJsonDnaHitOut(TConfig const& c, TStream& rcfile, std::vector<std::string> const& qn, std::vector<DnaHit> const& ht, std::vector<std::string> const& msg) {
    bool errors = false;
    rcfile << "{";
    // Errors
    rcfile << "\"errors\": [";
    for(uint32_t i = 0; i < msg.size(); ++i) {
      std::string msgtype = "warning";
      if (boost::starts_with(msg[i], "Error")) {	
	errors = true;
	msgtype = "error";
      }
      nlohmann::json err;
      err["type"] = msgtype;
      err["title"] = msg[i];
      if (i>0) rcfile << ',';
      rcfile << err.dump();
    }
    rcfile << "]";
    if (!errors) {
      // Meta information
      rcfile << ",\"meta\":";
      nlohmann::json meta;
      meta["version"] = diceyVersionNumber;
      meta["subcommand"] = "hunt";
      meta["distance"] = c.distance;
      meta["sequence"] = c.sequence;
      meta["genome"] = c.genome.string();
      meta["outfile"] = c.outfile.string();
      meta["maxmatches"] = c.max_locations;
      meta["hamming"] = (!c.indel);
      meta["forwardonly"] = (!c.reverse);
      rcfile << meta.dump() << ',';


      uint32_t oldchr = 999999;
      uint32_t oldstart = 0;
      rcfile << "\"data\":[";      
      for(uint32_t i = 0; i < ht.size(); ++i) {
	if ((oldchr != ht[i].chr) || (oldstart != ht[i].start)) {
	  if (i>0) rcfile << ',';
	  nlohmann::json j;
	  j["distance"] = std::abs(ht[i].score);
	  j["chr"] = qn[ht[i].chr];
	  j["start"] = ht[i].start;
	  j["end"] = ht[i].start + _nucleotideLength(ht[i].refalign) - 1;
	  j["strand"] = std::string(1, ht[i].strand);
	  j["refalign"] = ht[i].refalign;
	  j["queryalign"] = ht[i].queryalign;
	  rcfile << j.dump();
	}
	oldchr = ht[i].chr;
	oldstart = ht[i].start;
      }
      rcfile << ']';
    }
    rcfile << '}' << std::endl;
  }
  
  template<typename TConfig>
  inline void
  jsonDnaHitOut(TConfig const& c, std::vector<std::string> const& qn, std::vector<DnaHit> const& ht, std::vector<std::string> const& msg) {
    if (c.hasOutfile) {
      // Output file
      boost::iostreams::filtering_ostream rcfile;
      rcfile.push(boost::iostreams::gzip_compressor());
      rcfile.push(boost::iostreams::file_sink(c.outfile.c_str(), std::ios_base::out | std::ios_base::binary));
      writeJsonDnaHitOut(c, rcfile, qn, ht, msg);
      rcfile.pop();
    } else {
      writeJsonDnaHitOut(c, std::cout, qn, ht, msg);
    }
  }
  
  int hunter(int argc, char** argv) {
    HunterConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "gzipped output file")
      ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(1000), "max. number of matches")
      ("maxNeighborhood,x", boost::program_options::value<uint32_t>(&c.maxNeighborhood)->default_value(10000), "max. neighborhood size")
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
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] -g Danio_rerio.fa.gz CATTACTAACATCAGT" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // DNA Hits
    std::vector<DnaHit> ht;
    std::vector<std::string> msg;
    std::vector<uint32_t> seqlen;
    std::vector<std::string> seqname;
    
    // Cmd switches
    if (!vm.count("hamming")) c.indel = true;
    else c.indel = false;
    if (!vm.count("forward")) c.reverse = true;
    else c.reverse = false;
    if (vm.count("outfile")) c.hasOutfile = true;
    else c.hasOutfile = false;

    // Check sequence length
    if (c.sequence.size() < 10) {
      msg.push_back("Error: Input sequence is shorter than 10 nucleotides!");
      jsonDnaHitOut(c, seqname, ht, msg);
      return 1;
    }

    // Upper case
    c.sequence = boost::to_upper_copy(c.sequence);
    c.sequence = replaceNonDna(c.sequence, msg);
    std::string revSequence = c.sequence;
    reverseComplement(revSequence);
    
    // Check distance
    if (c.distance >= c.sequence.size()) {
      c.distance = c.sequence.size() - 1;
      msg.push_back("Warning: Distance was adjusted to sequence length!");
    }

    // Set prefix and suffix based on edit distance
    c.pre_context = 0;
    c.post_context = 0;
    if (c.indel) {
      c.pre_context += c.distance;
      c.post_context += c.distance;
    }
    
    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      msg.push_back("Error: Genome does not exist!");
      jsonDnaHitOut(c, seqname, ht, msg);
      return 1;
    }
    
    // Parse chromosome lengths
    uint32_t nseq = getSeqLenName(c, seqlen, seqname);
    if (!nseq) {
      msg.push_back("Error: Could not retrieve sequence lengths!");
      jsonDnaHitOut(c, seqname, ht, msg);
      return 1;
    }

    // Reference index
    csa_wt<> fm_index;  
    boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
    std::string index_file = op.string() + ".fm9";
    if (!load_from_file(fm_index, index_file)) {
      msg.push_back("Error: FM-Index cannot be loaded!");
      jsonDnaHitOut(c, seqname, ht, msg);
      return 1;
    }

    // Define alphabet
    typedef std::set<char> TAlphabet;
    char tmp[] = {'A', 'C', 'G', 'T'};
    TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));

    // Generate neighbors
    typedef std::set<std::string> TStringSet;
    typedef std::vector<TStringSet> TFwdRevSearchSets;
    TFwdRevSearchSets fwrv(2, TStringSet());
    neighbors(c.sequence, alphabet, c.distance, c.indel, c.maxNeighborhood, fwrv[0]);
    // Debug
    //for(TStringSet::iterator it = fwrv[0].begin(); it != fwrv[0].end(); ++it) std::cerr << *it << std::endl;

    // Reverse complement set?
    if (c.reverse) neighbors(revSequence, alphabet, c.distance, c.indel, c.maxNeighborhood, fwrv[1]);

    // Check neighborhood size
    if ((fwrv[0].size() >= c.maxNeighborhood) || (fwrv[1].size() >= c.maxNeighborhood)) {
      std::string m = "Warning: Neighborhood size exceeds " + boost::lexical_cast<std::string>(c.maxNeighborhood) + " candidates. Only first " + boost::lexical_cast<std::string>(c.maxNeighborhood) + " neighbors are searched, results are likely incomplete!";
      msg.push_back(m);
    }
    
    // Search
    uint32_t hits = 0;
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
	    std::size_t pre_extract = c.pre_context;
	    std::size_t post_extract = c.post_context;
	    if (pre_extract > locations[i]) {
	      pre_extract = locations[i];
	    }
	    if (locations[i]+m+post_extract > fm_index.size()) {
	      post_extract = fm_index.size() - locations[i] - m;
	    }
	    auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
	    std::string pre = s.substr(0, pre_extract);
	    s = s.substr(pre_extract);
	    if (pre.find_last_of('\n') != std::string::npos) {
	      pre = pre.substr(pre.find_last_of('\n')+1);
	    }
	    std::string post = s.substr(m);
	    post = post.substr(0, post.find_first_of('\n'));
	    
	    // Genomic sequence
	    std::string genomicseq = pre + s.substr(0, m) + post;
	    if (pre.size() < chrpos) chrpos -= pre.size();
	    DnaScore<int32_t> sc(0, -1, -1, -1);
	    typedef boost::multi_array<char, 2> TAlign;
	    AlignConfig<false, true> global;
	    if (fwrvidx == 0) {
	      if (c.indel) {
		TAlign align;
		int32_t score = needle(genomicseq, c.sequence , align, global, sc);
		std::string refalign = "";
		std::string queryalign = "";
		bool leadGap = true;
		for(uint32_t j = 0; (j < (align.shape()[1] - _trailGap(align))); ++j) {
		  if (align[1][j] != '-') leadGap = false;
		  if (!leadGap) {
		    refalign += align[0][j];
		    queryalign += align[1][j];
		  } else {
		    ++chrpos;
		  }
		}
		ht.push_back(DnaHit(score, refIndex, chrpos+1, '+', refalign, queryalign));
	      } else {
		int32_t score = needleScore(genomicseq, c.sequence, sc);
		ht.push_back(DnaHit(score, refIndex, chrpos+1, '+', genomicseq, c.sequence));
	      }
	    } else {
	      if (c.indel) {
		TAlign align;
		int32_t score = needle(genomicseq, revSequence, align, global, sc);
		std::string refalign = "";
		std::string queryalign = "";
		bool leadGap = true;
		for(uint32_t j = 0; (j < (align.shape()[1] - _trailGap(align))); ++j) {
		  if (align[1][j] != '-') leadGap = false;
		  if (!leadGap) {
		    refalign += align[0][j];
		    queryalign += align[1][j];
		  } else {
		    ++chrpos;
		  }
		}
		ht.push_back(DnaHit(score, refIndex, chrpos+1, '-', refalign, queryalign));
	      } else {
		int32_t score = needleScore(genomicseq, revSequence, sc);
		ht.push_back(DnaHit(score, refIndex, chrpos+1, '-', genomicseq, revSequence));
	      }
	    }
	    ++hits;
	  }
	}
      }
    }
    if (hits >= c.max_locations) {
      std::string m = "Warning: More than " + boost::lexical_cast<std::string>(c.max_locations) + " matches found. Only first " + boost::lexical_cast<std::string>(c.max_locations) + " matches are reported, results are likely incomplete!";
      msg.push_back(m);
    }
    
    // Sort
    std::sort(ht.begin(), ht.end(), SortDnaHit<DnaHit>());

    // Output
    jsonDnaHitOut(c, seqname, ht, msg);
    
    return 0;
  }

}

#endif
