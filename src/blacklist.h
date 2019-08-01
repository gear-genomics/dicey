#ifndef BLACKLIST_H
#define BLACKLIST_H

#include <boost/icl/split_interval_map.hpp>
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

  struct BlacklistConfig {
    boost::filesystem::path map;
    boost::filesystem::path blacklist;
    boost::filesystem::path outfile;
  };


  int blacklist(int argc, char** argv) {
    BlacklistConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("blacklist,b", boost::program_options::value<boost::filesystem::path>(&c.blacklist)->default_value("blacklist.bed"), "blacklist in BED format")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("map.fa.gz"), "gzipped output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.map), "mappability map")
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
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] -b blacklist.bed Danio_rerio.fa.gz" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Check mappability map
    if (!(boost::filesystem::exists(c.map) && boost::filesystem::is_regular_file(c.map) && boost::filesystem::file_size(c.map))) {
      std::cerr << "Error: Mappability map does not exist!" << std::endl;
      return 1;
    }

    // Check blacklist
    if (!(boost::filesystem::exists(c.blacklist) && boost::filesystem::is_regular_file(c.blacklist) && boost::filesystem::file_size(c.blacklist))) {
      std::cerr << "Error: Blacklist does not exist!" << std::endl;
      return 1;
    }

    // Outfile
    boost::iostreams::filtering_ostream of;
    of.push(boost::iostreams::gzip_compressor());
    of.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Get chromosomes
    faidx_t* fai = fai_load(c.map.string().c_str());
    uint32_t nchr = faidx_nseq(fai);
    std::vector<std::string> chrName(nchr);
    for(uint32_t refIndex = 0; refIndex < nchr; ++refIndex) chrName[refIndex] = std::string(faidx_iseq(fai, refIndex));

    // Get blacklist regions
    typedef boost::icl::interval_set<uint32_t> TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome blacklistRegions;
    if (!_parseBedIntervals(c.blacklist.string(), chrName, blacklistRegions)) {
      std::cerr << "Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
      return 1;
    }
   
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Blacklist mappability map" << std::endl;
    boost::progress_display show_progress( nchr );
    for(uint32_t refIndex = 0; refIndex < nchr; ++refIndex) {
      ++show_progress;

      // Load chromosome
      std::string seqname(faidx_iseq(fai, refIndex));
      int32_t sql = faidx_seq_len(fai, seqname.c_str());
      int32_t seqlen = -1;
      char* seq = faidx_fetch_seq(fai, seqname.c_str(), 0, sql, &seqlen);

      // Blacklist
      for(typename TChrIntervals::iterator it = blacklistRegions[refIndex].begin(); it != blacklistRegions[refIndex].end(); ++it) {
	if ((it->lower() < it->upper()) && ((int32_t) it->upper() < sql)) {
	  for(uint32_t k = it->lower(); k < it->upper(); ++k) {
	    seq[k] = 'N';
	  }
	}
      }

      // Output
      of << ">" << seqname << std::endl;
      for(int32_t i = 0; i<sql; ++i) of << seq[i];
      of << std::endl;

      // Clean-up	
      if (seq != NULL) free(seq);
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
