#ifndef MAPBAM_H
#define MAPBAM_H

#include <sdsl/suffix_arrays.hpp>

#include <fstream>
#include <iomanip>

#include <boost/unordered_map.hpp>
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
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "neighbors.h"
#include "util.h"

using namespace sdsl;

namespace dicey
{

  struct MapBamConfig {
    bool hasChr;
    uint16_t minQual;
    int32_t isize;
    std::string chrom;
    boost::filesystem::path bamFile;
    boost::filesystem::path outfile;
  };


  int mapbam(int argc, char** argv) {
    MapBamConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("chromosome,c", boost::program_options::value<std::string>(&c.chrom), "chromosome name to process")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("map.fa.gz"), "gzipped output file")
      ("insertsize,s", boost::program_options::value<int32_t>(&c.isize)->default_value(501), "insert size")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "BAM file")
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
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] chopped.bam" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Chromosome
    if (vm.count("chromosome")) c.hasChr = true;
    else c.hasChr = false;
    
    // Half-window
    int32_t halfwin = (int32_t) (c.isize / 2);
    c.isize = 2 * halfwin + 1;

    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Outfile
    boost::iostreams::filtering_ostream of;
    of.push(boost::iostreams::gzip_compressor());
    of.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Mappability" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    // Iterate chromosomes
    for(uint32_t refIndex = 0; refIndex < (uint32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if ((c.hasChr) && (std::string(hdr->target_name[refIndex]) != c.chrom)) continue;
      if (chrNoData(c, refIndex, idx)) continue;

      // Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov(hdr->target_len[refIndex], 0);
      TCoverage ucov(hdr->target_len[refIndex], 0);
      
      // Mate map
      typedef boost::unordered_map<std::size_t, bool> TMateMap;
      TMateMap mateMap;
      
      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((rec->core.flag & BAM_FPAIRED) && ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid))) continue;

	int32_t midPoint = rec->core.pos + halfAlignmentLength(rec);
	int32_t isize = 0;
	if (rec->core.flag & BAM_FPAIRED) {
	  // Clean-up the read store for identical alignment positions
	  if (rec->core.pos > lastAlignedPos) {
	    lastAlignedPosReads.clear();
	    lastAlignedPos = rec->core.pos;
	  }

	  if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	    // First read
	    lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	    std::size_t hv = hash_pair(rec);
	    mateMap[hv] = true;
	    continue;
	  } else {
	    // Second read
	    std::size_t hv = hash_pair_mate(rec);
	    if ((mateMap.find(hv) == mateMap.end()) || (!mateMap[hv])) continue; // Mate discarded
	    mateMap[hv] = false;
	  }
	  isize = (rec->core.pos + alignmentLength(rec)) - rec->core.mpos;
	  midPoint = rec->core.mpos + (int32_t) (isize/2);
	}

	// Count fragment
	if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (cov[midPoint] < maxCoverage - 1)) ++cov[midPoint];
	if ((rec->core.qual >= c.minQual) && (isize == c.isize)) {
	  if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex]) && (ucov[midPoint] < maxCoverage - 1)) ++ucov[midPoint];
	}
      }
      bam_destroy1(rec);
      hts_itr_destroy(iter);

      // Fill map
      of << ">" << hdr->target_name[refIndex] << std::endl;
      for(uint32_t pos = 0; pos < hdr->target_len[refIndex]; ++pos) {
	if ((cov[pos] == 0) || (cov[pos] > 1)) of << 'N';
	else {
	  if (cov[pos] != ucov[pos]) of << 'A';
	  else of << 'C';
	}
      }
      of << std::endl;
    }

    // Clean-up
    of.pop();

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
