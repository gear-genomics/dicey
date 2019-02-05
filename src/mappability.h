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

#ifndef MAPPABILITY_H
#define MAPPABILITY_H

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

#ifdef OPENMP
#include <omp.h>
#endif

#include <htslib/faidx.h>

#include "neighbors.h"
#include "util.h"

using namespace sdsl;

namespace dicey
{

  struct MapConfig {
    uint32_t maxEditDistance;
    uint32_t maxNeighborhood;
    int32_t readlength;
    int32_t chrom;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
  };


  int mappability(int argc, char** argv) {
    MapConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("map.fa.gz"), "gzipped output file")
      ("readlength,r", boost::program_options::value<int32_t>(&c.readlength)->default_value(51), "read length")
      ("chromosome,c", boost::program_options::value<int32_t>(&c.chrom)->default_value(-1), "chromosome index (-1: all)")
      ("maxEditDistance,e", boost::program_options::value<uint32_t>(&c.maxEditDistance)->default_value(1), "max. edit distance")
      ("maxNeighborhood,x", boost::program_options::value<uint32_t>(&c.maxNeighborhood)->default_value(10000), "max. neighborhood size")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.genome), "indexed genome")
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
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] Danio_rerio.fa.gz" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Half-window
    int32_t halfwin = (int32_t) (c.readlength / 2);
    c.readlength = 2 * halfwin + 1;

    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Error: Genome does not exist!" << std::endl;
      return 1;
    }
    
    // Reference index
    csa_wt<> fm_index;  
    boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
    std::string index_file = op.string() + ".fm9";
    if (!load_from_file(fm_index, index_file)) {
      std::cerr << "Error: FM-Index cannot be loaded!" << std::endl;
      return 1;
    }

    // Define alphabet
    typedef std::set<char> TAlphabet;
    char tmp[] = {'A', 'C', 'G', 'T'};
    TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
    
    // Iterate chromosomes
    faidx_t* fai = fai_load(c.genome.string().c_str());
    uint32_t nchr = faidx_nseq(fai);
    if ((c.chrom != -1) && (c.chrom >= (int32_t) nchr)) {
      std::cerr << "Only " << nchr << " chromosomes in reference genome!" << std::endl;
      fai_destroy(fai);
      return 1;
    }

    // Outfile
    boost::iostreams::filtering_ostream of;
    of.push(boost::iostreams::gzip_compressor());
    of.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Mappability" << std::endl;
    boost::progress_display show_progress( nchr );
    
    int64_t idxpos = 0;
    for(uint32_t refIndex = 0; refIndex < nchr; ++refIndex) {
      ++show_progress;
      if ((c.chrom != -1) && ((int32_t) refIndex != c.chrom)) continue;
      int32_t seqlen = -1;
      std::string seqname(faidx_iseq(fai, refIndex));
      char* seq = faidx_fetch_seq(fai, seqname.c_str(), 0, faidx_seq_len(fai, seqname.c_str()), &seqlen);
      std::vector<char> oseq(seqlen, 'N');
#pragma omp parallel for default(shared)
      for(int32_t pos = 0; pos < (seqlen - c.readlength + 1); ++pos) {
	int64_t ipos = idxpos + pos;
	std::string sequence = boost::to_upper_copy(std::string(seq + pos, seq + pos + c.readlength));
	bool nContent = false;
	for(uint32_t i = 0; i < sequence.size(); ++i) {
	  if (sequence[i] == 'N') {
	    nContent = true;
	    break;
	  }
	}
	if (!nContent) {
	    // Reverse complement
	  std::string revSequence = sequence;
	  reverseComplement(revSequence);
	  
	  // Find max. unique edit distance
	  int32_t maxED = -1;
	  bool unique = true;
	  for(uint32_t ed = 0; ((ed <= c.maxEditDistance) && (unique)); ++ed) {
	    // Neighbors
	    typedef std::set<std::string> TStringSet;
	    typedef std::vector<TStringSet> TFwdRevSearchSets;
	    TFwdRevSearchSets fwrv(2, TStringSet());
	    neighbors(sequence, alphabet, ed, true, c.maxNeighborhood, fwrv[0]);
	    neighbors(revSequence, alphabet, ed, true, c.maxNeighborhood, fwrv[1]);
	    
	    // Check neighborhood size
	    if ((fwrv[0].size() >= c.maxNeighborhood) || (fwrv[1].size() >= c.maxNeighborhood)) {
	      std::cerr << "Warning: Neighborhood size exceeds " + boost::lexical_cast<std::string>(c.maxNeighborhood) + " candidates. Only first " + boost::lexical_cast<std::string>(c.maxNeighborhood) + " neighbors are searched, results are likely incomplete!" << std::endl;
	    }
	    
	    // Search
	    for(uint32_t fwrvidx = 0; ((fwrvidx < fwrv.size()) && (unique)); ++fwrvidx) {
	      for(typename TStringSet::const_iterator it = fwrv[fwrvidx].begin(); it != fwrv[fwrvidx].end(); ++it) {
		std::string query = *it;
		std::size_t m = query.size();
		std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
		if (occs > 1) {
		  unique = false;
		  break;
		  } else {
		  if ((occs) && (ed)) {
		    // Inexact matches, check position is correct
		    auto locations = locate(fm_index, query.begin(), query.begin() + m);
		    if (std::abs((int64_t) locations[0] - ipos) > ed) {
		      unique = false;
		      break;
		    }
		  }
		  }
	      }
	    }
	    
	    // Accept edit distance?
	    if (unique) {
	      maxED = ed;
	    }
	  }
	  if (maxED == 0) oseq[pos+halfwin]='A';
	  else if (maxED == 1) oseq[pos+halfwin]='C';
	  else if (maxED == 2) oseq[pos+halfwin]='G';
	  else if (maxED > 2) oseq[pos+halfwin]='T';
	}
      }
      of << ">" << seqname << std::endl;
      for(uint32_t i = 0; i<oseq.size(); ++i) of << oseq[i];
      of << std::endl;
      if (seq != NULL) free(seq);
      idxpos += seqlen + 1;
    }

    // Clean-up
    fai_destroy(fai);
    of.pop();

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
