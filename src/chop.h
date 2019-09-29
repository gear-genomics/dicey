#ifndef CHOP_H
#define CHOP_H

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

  struct ChopConfig {
    bool se;
    bool sf;
    bool revcomp;
    bool hashing;
    uint32_t readlength;
    uint32_t isize;
    boost::filesystem::path genome;
    std::string fq1;
    std::string fq2;
  };


  template<typename TConfig, typename TStream>
  inline void
  _processFasta(TConfig const& c, TStream& of, std::string const& seq) {
    for(uint32_t pos = 0; ((pos + c.readlength) <= seq.size()); ++pos) {
      std::string kmer = boost::to_upper_copy(std::string(&seq[0] + pos, &seq[0] + pos + c.readlength));
      if (nContent(kmer)) continue;
      unsigned h1 = hash_string(kmer.c_str());
      reverseComplement(kmer);
      unsigned h2 = hash_string(kmer.c_str());
      if (h1 < h2) of << h1 << '\t' << h2 << std::endl;
      else of << h2 << '\t' << h1 << std::endl;
    }
  }
  
  int chop(int argc, char** argv) {
    ChopConfig c;
    std::string qual("AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJAJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJ7FJJJFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("fq1,f", boost::program_options::value<std::string>(&c.fq1)->default_value("read1"), "read1 output prefix")
      ("fq2,g", boost::program_options::value<std::string>(&c.fq2)->default_value("read2"), "read2 output prefix")
      ("length,l", boost::program_options::value<uint32_t>(&c.readlength)->default_value(101), "read length")
      ("insertsize,i", boost::program_options::value<uint32_t>(&c.isize)->default_value(501), "insert size")
      ("se,s", "generate single-end data")
      ("chromosome,c", "generate reads by chromosome")
      ("revcomp,r", "reverse complement all reads")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.genome), "indexed genome")
      ("hashing,a", "output hashed reads")
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

    // Single-end?
    if (vm.count("se")) c.se = true;
    else c.se = false;

    // Single output file?
    if (vm.count("chromosome")) c.sf = false;
    else c.sf = true;

    // Reverse complement
    if (vm.count("revcomp")) c.revcomp = true;
    else c.revcomp = false;

    // Hashing
    if (vm.count("hashing")) {
      // Kmer mode
      c.hashing = true;
      c.se = true; // Always single-end, fwd and rev
    } else c.hashing = false;

    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Error: Genome does not exist!" << std::endl;
      return 1;
    }

    // Single-end or paired-end?
    if (c.se) {
      // Read-length
      if (c.readlength >= qual.size()) c.readlength = qual.size() - 1;
      std::string readQual = qual.substr(0, c.readlength);

      // Hashing mode?
      if (c.hashing) {
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Kmer hashes" << std::endl;
	
	// Output file
	boost::iostreams::filtering_ostream of1;
	std::string read1fq = c.fq1 + ".hashes.gz";
	of1.push(boost::iostreams::gzip_compressor());
	of1.push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));

	// Fasta input (gzipped)
	std::string faname = "";
	std::string tmpfasta = "";
	std::ifstream fafile(c.genome.string().c_str(), std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	dataIn.push(boost::iostreams::gzip_decompressor());
	dataIn.push(fafile);
	std::istream instream(&dataIn);
	std::string line;
	while(std::getline(instream, line)) {
	  if (!line.empty()) {
	    if (line[0] == '>') {
	      if (!tmpfasta.empty()) {
		_processFasta(c, of1, tmpfasta);
		tmpfasta.clear();
	      }
	      if (line.at(line.length() - 1) == '\r' ){
		faname = line.substr(1, line.length() - 2);
	      } else {
		faname = line.substr(1);
	      }
	    } else {
	      if (line.at(line.length() - 1) == '\r' ){
		tmpfasta += boost::to_upper_copy(line.substr(0, line.length() - 1));
	      } else {
		tmpfasta += boost::to_upper_copy(line);
	      }
	    }
	  }
	}
	if (!tmpfasta.empty()) _processFasta(c, of1, tmpfasta);
	dataIn.pop();
	dataIn.pop();
	fafile.close();

	of1.pop();
	of1.pop();
      } else {
	// Iterate chromosomes
	faidx_t* fai = fai_load(c.genome.string().c_str());
	uint32_t nchr = faidx_nseq(fai);
      
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Chop reference" << std::endl;
	boost::progress_display show_progress( nchr );

	uint32_t index = 0;
	boost::iostreams::filtering_ostream of1;
	if (c.sf) {
	  std::string read1fq = c.fq1 + ".fq.gz";
	  of1.push(boost::iostreams::gzip_compressor());
	  of1.push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));
	}
	for(uint32_t refIndex = 0; refIndex < nchr; ++refIndex) {
	  ++show_progress;

	  // Load chromosome
	  std::string seqname(faidx_iseq(fai, refIndex));
	  int32_t sql = faidx_seq_len(fai, seqname.c_str());
	  int32_t seqlen = -1;
	  char* seq = faidx_fetch_seq(fai, seqname.c_str(), 0, sql, &seqlen);
	
	  // Generate paired-end reads
	  if ((int32_t) c.readlength <= sql) {
	  
	    // FQ1
	    if (!c.sf) {
	      std::string read1fq = c.fq1 + "." + seqname + ".fq.gz";
	      of1.push(boost::iostreams::gzip_compressor());
	      of1.push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));
	    }
	  
	    // Iterate chr
	    for(int32_t pos = 0; ((pos + (int32_t) c.readlength) <= sql); ++pos, ++index) {
	      std::string read1 = boost::to_upper_copy(std::string(seq + pos, seq + pos + c.readlength));
	      if (nContent(read1)) continue;
	      if (c.revcomp) reverseComplement(read1);
	      of1 << "@Frag" << index << "_" << seqname << "_" << pos << " 1:N:0:0" << std::endl;
	      of1 << read1 << std::endl;
	      of1 << "+" << std::endl;
	      of1 << readQual << std::endl;
	    }
	    if (!c.sf) {
	      of1.pop();
	      of1.pop();
	    }
	  }
	  // Clean-up	
	  if (seq != NULL) free(seq);
	}
	if (c.sf) {
	  of1.pop();
	  of1.pop();
	}

	// Clean-up
	fai_destroy(fai);
      }
      
    } else {
      // Half-window
      int32_t halfwin = (int32_t) (c.isize / 2);
      c.isize = 2 * halfwin + 1;

      // Read-length < isize
      if (c.readlength > c.isize) c.readlength = c.isize - 1;
      if (c.readlength >= qual.size()) c.readlength = qual.size() - 1;
      std::string readQual = qual.substr(0, c.readlength);

      // Iterate chromosomes
      faidx_t* fai = fai_load(c.genome.string().c_str());
      uint32_t nchr = faidx_nseq(fai);
      
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Chop reference" << std::endl;
      boost::progress_display show_progress( nchr );

      uint32_t index = 0;
      boost::iostreams::filtering_ostream of1;
      if (c.sf) {
	std::string read1fq = c.fq1 + ".fq.gz";
	of1.push(boost::iostreams::gzip_compressor());
	of1.push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));
      }
      boost::iostreams::filtering_ostream of2;
      if (c.sf) {
	std::string read2fq = c.fq2 + ".fq.gz";
	of2.push(boost::iostreams::gzip_compressor());
	of2.push(boost::iostreams::file_sink(read2fq.c_str(), std::ios_base::out | std::ios_base::binary));
      }
      for(uint32_t refIndex = 0; refIndex < nchr; ++refIndex) {
	++show_progress;

	// Load chromosome
	std::string seqname(faidx_iseq(fai, refIndex));
	int32_t sql = faidx_seq_len(fai, seqname.c_str());
	int32_t seqlen = -1;
	char* seq = faidx_fetch_seq(fai, seqname.c_str(), 0, sql, &seqlen);
	
	// Generate paired-end reads
	if (halfwin < sql) {
	  
	  // FQ1
	  if (!c.sf) {
	    std::string read1fq = c.fq1 + "." + seqname + ".fq.gz";
	    of1.push(boost::iostreams::gzip_compressor());
	    of1.push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));
	  }
	  
	  // FQ2
	  if (!c.sf) {
	    std::string read2fq = c.fq2 + "." + seqname + ".fq.gz";
	    of2.push(boost::iostreams::gzip_compressor());
	    of2.push(boost::iostreams::file_sink(read2fq.c_str(), std::ios_base::out | std::ios_base::binary));
	  }
	  
	  // Iterate chr
	  for(int32_t pos = halfwin; pos < sql - halfwin; ++pos, ++index) {
	    std::string read1 = boost::to_upper_copy(std::string(seq + pos - halfwin, seq + pos - halfwin + c.readlength));
	    if (nContent(read1)) continue;
	    if (c.revcomp) reverseComplement(read1);
	    std::string read2 = boost::to_upper_copy(std::string(seq + pos + halfwin - c.readlength + 1, seq + pos + halfwin + 1));
	    if (nContent(read2)) continue;
	    if (!c.revcomp) reverseComplement(read2);
	    of1 << "@Frag" << index << "_" << seqname << "_" << pos << " 1:N:0:0" << std::endl;
	    of1 << read1 << std::endl;
	    of1 << "+" << std::endl;
	    of1 << readQual << std::endl;
	    of2 << "@Frag" << index << "_" << seqname << "_" << pos << " 2:N:0:0" << std::endl;
	    of2 << read2 << std::endl;
	    of2 << "+" << std::endl;
	    of2 << readQual << std::endl;
	  }
	  if (!c.sf) {
	    of1.pop();
	    of1.pop();
	    of2.pop();
	    of2.pop();
	  }
	}
	// Clean-up	
	if (seq != NULL) free(seq);
      }
      if (c.sf) {
	of1.pop();
	of1.pop();
	of2.pop();
	of2.pop();
      }

      // Clean-up
      fai_destroy(fai);
    }
      
    // Done
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
