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

#include "neighbors.h"
#include "util.h"

using namespace sdsl;

namespace dicey
{

  struct HunterConfig {
    bool indel;
    bool reverse;
    uint32_t distance;
    std::size_t pre_context;
    std::size_t post_context;
    std::size_t max_locations;
    std::string sequence;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
  };

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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] AGGTACTAACG" << std::endl;
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

    // Check distance
    if (c.distance >= c.sequence.size()) c.distance = c.sequence.size() - 1;

    // Set prefix and suffix based on edit distance
    c.pre_context = 0;
    c.post_context = 0;
    if (c.indel) {
      c.pre_context += c.distance;
      c.post_context += c.distance;
    }

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
    //for(TStringSet::iterator it = fwdset.begin(); it != fwdset.end(); ++it) std::cerr << *it << std::endl;

    // Reverse complement set?
    if (c.reverse) {
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Generate reverse sequence neighborhood " << std::flush;
      std::string rev = c.sequence;
      reverseComplement(rev);
      neighbors(rev, alphabet, c.distance, c.indel, fwrv[1]);
      std::cout << "(#n=" << fwrv[1].size() << ")" << std::endl;
    }

    // Search
    uint32_t hits = 0;
    for(uint32_t i = 0; i < fwrv.size(); ++i) {
      for(typename TStringSet::const_iterator it = fwrv[i].begin(); ((it != fwrv[i].end()) && (hits < c.max_locations)); ++it) {
	std::string query = *it;
	std::size_t m = query.size();
	std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
	if (occs > 0) {
	  auto locations = locate(fm_index, query.begin(), query.begin() + m);
	  std::sort(locations.begin(), locations.end());
	  std::size_t pre_extract = c.pre_context;
	  std::size_t post_extract = c.post_context;
	  for(std::size_t i = 0; ((i < std::min(occs, c.max_locations)) && (hits < c.max_locations)); ++i) {
	    ++hits;
	    int64_t bestPos = locations[i];
	    int64_t cumsum = 0;
	    uint32_t refIndex = 0;
	    for(; bestPos >= cumsum + seqlen[refIndex]; ++refIndex) cumsum += seqlen[refIndex];
	    uint32_t chrpos = bestPos - cumsum;
	    if (pre_extract > locations[i]) {
	      pre_extract = locations[i];
	    }
	    if (locations[i]+m+post_extract > fm_index.size()) {
	      post_extract = fm_index.size() - locations[i] - m;
	    }
	    auto s = extract(fm_index, locations[i]-pre_extract, locations[i]+m+post_extract-1);
	    std::string pre = s.substr(0, pre_extract);
	    s = s.substr(pre_extract);
	    if (pre.find_last_of('\n') != std::string::npos) {
	      pre = pre.substr(pre.find_last_of('\n')+1);
	    }
	    std::string post = s.substr(m);
	    post = post.substr(0, post.find_first_of('\n'));
	    std::string genomicseq = pre + s.substr(0, m) + post;
	    std::cout << seqname[refIndex] << ":" << chrpos+1 << "-" << chrpos+query.size() << std::endl;
	    std::cout << genomicseq << std::endl;
	    std::cout << c.sequence << std::endl;
	    std::cout << "--" << std::endl;
	  }
	}
      }
    }
    
    return 0;
  }

}

#endif
