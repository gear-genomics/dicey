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
    bool dumphash;
    uint32_t readlength;
    uint32_t isize;
    uint32_t nTmpFile;
    boost::filesystem::path genome;
    std::string fq1;
    std::string fq2;
  };


  template<typename TConfig, typename TStream>
  inline void
  _processFasta(TConfig const& c, TStream& of, std::string const& seq, uint64_t& kmerCount) {
    std::string rcseq(seq);
    reverseComplement(rcseq);
    uint32_t seqlen = seq.size();
    for(uint32_t pos = 0; ((pos + c.readlength) <= seqlen); ++pos) {
      if (nContent(seq.substr(pos, c.readlength))) continue;
      unsigned h1 = hash_string(seq.substr(pos, c.readlength).c_str());
      unsigned h2 = hash_string(rcseq.substr(seqlen - c.readlength - pos, c.readlength).c_str());
      if (h1 < h2) {
	//std::cerr << h1 << '\t' << h2 << std::endl;
	of << h1 << '\t' << h2 << std::endl;
	if (c.dumphash) {
	  std::cerr << h1 << '\t' << h2 << '\t' << seq.substr(pos, c.readlength) << '\t' << rcseq.substr(seqlen - c.readlength - pos, c.readlength) << std::endl;
	}
      } else {
	//std::cerr << h2 << '\t' << h1 << std::endl;
	of << h2 << '\t' << h1 << std::endl;
	if (c.dumphash) {
	  std::cerr << h2 << '\t' << h1 << '\t' << rcseq.substr(seqlen - c.readlength - pos, c.readlength) << '\t' << seq.substr(pos, c.readlength) << std::endl;
	}
      }
      ++kmerCount;
    }
  }
  
  int chop(int argc, char** argv) {
    ChopConfig c;
    c.nTmpFile = 100;
    
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
      ("dumphash,d", "dump hashed k-mers")
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
      if (vm.count("dumphash")) c.dumphash = true;
      else c.dumphash = false;
    } else c.hashing = false;

    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Error: Genome does not exist!" << std::endl;
      return 1;
    }

    // Single-end or paired-end?
    if (c.se) {
      // Read-length
      while (c.readlength >= qual.size()) qual.push_back('A');
      std::string readQual = qual.substr(0, c.readlength);

      // Hashing mode?
      if (c.hashing) {
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Kmer hashes" << std::endl;
	
	// Output file
	uint64_t lastKmerCount = 0;
	uint64_t kmerCount = 0;
	uint32_t fileIdx = 0;
	std::vector<boost::iostreams::filtering_ostream> ofAll(c.nTmpFile);
	for(uint32_t i = 0; i < c.nTmpFile; ++i) {
	  std::string read1fq = c.fq1 + "." + boost::lexical_cast<std::string>(i) + ".hashes.gz";
	  ofAll[i].push(boost::iostreams::gzip_compressor());
	  ofAll[i].push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));
	}

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
	  // Swap every 1M hashes chunk file
	  if (kmerCount - lastKmerCount > 1000000) {
	    lastKmerCount = kmerCount;
	    ++fileIdx;
	    if (fileIdx >= c.nTmpFile) fileIdx = 0;
	  }
	  //if (kmerCount > 1000030) break;
	  if (!line.empty()) {
	    if (line[0] == '>') {
	      if (!tmpfasta.empty()) {
		_processFasta(c, ofAll[fileIdx], tmpfasta, kmerCount);
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
	if (!tmpfasta.empty()) _processFasta(c, ofAll[fileIdx], tmpfasta, kmerCount);
	dataIn.pop();
	dataIn.pop();
	fafile.close();
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Hashed " << kmerCount << " k-mers." << std::endl;	

	for(uint32_t i= 0; i < c.nTmpFile; ++i) {
	  ofAll[i].pop();
	  ofAll[i].pop();
	}

	// Sort chunks
	for(uint32_t i = 0; i < c.nTmpFile; ++i) {
	  std::string read1fq = c.fq1 + "." + boost::lexical_cast<std::string>(i) + ".hashes.gz";
	  std::ifstream kmerchunk(read1fq.c_str());
	  boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	  dataIn.push(boost::iostreams::gzip_decompressor());
	  dataIn.push(kmerchunk);
	  std::istream instream(&dataIn);
	  std::string line;
	  std::vector<std::pair<unsigned, unsigned> > hashVec;
	  while(std::getline(instream, line)) {
	    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	    boost::char_separator<char> sep("\t");
	    Tokenizer tokens(line, sep);
	    Tokenizer::iterator tokIter = tokens.begin();
	    unsigned h1 = boost::lexical_cast<unsigned>(*tokIter++);
	    unsigned h2 = boost::lexical_cast<unsigned>(*tokIter++);
	    hashVec.push_back(std::make_pair(h1, h2));
	  }
	  dataIn.pop();
	  dataIn.pop();
	  
	  // Delete input file
	  boost::filesystem::remove(read1fq);
	  
	  // Sort
	  std::sort(hashVec.begin(), hashVec.end());

	  // Sorted output file
	  boost::iostreams::filtering_ostream ofSort;
	  std::string read1sort = c.fq1 + "." + boost::lexical_cast<std::string>(i) + ".sort.gz";
	  ofSort.push(boost::iostreams::gzip_compressor());
	  ofSort.push(boost::iostreams::file_sink(read1sort.c_str(), std::ios_base::out | std::ios_base::binary));
	  for(uint32_t j = 0; j < hashVec.size(); ++j) {
	    ofSort << hashVec[j].first << '\t' << hashVec[j].second << std::endl;
	  }
	  ofSort.pop();
	  ofSort.pop();
	}
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Sorted " << c.nTmpFile << " chunks." << std::endl;	

	// Open final output file
	boost::iostreams::filtering_ostream dataOut;
	std::string read1fq = c.fq1 + ".hashes.gz";
	dataOut.push(boost::iostreams::gzip_compressor());
	dataOut.push(boost::iostreams::file_sink(read1fq.c_str(), std::ios_base::out | std::ios_base::binary));
	
	// Merge sorted chunks
	std::vector<std::ifstream> ifData(c.nTmpFile);
	std::vector<boost::iostreams::filtering_streambuf<boost::iostreams::input> > inData(c.nTmpFile);
	for(uint32_t i = 0; i < c.nTmpFile; ++i) {
	  std::string read1sort = c.fq1 + "." + boost::lexical_cast<std::string>(i) + ".sort.gz";
	  ifData[i].open(read1sort.c_str());
	  inData[i].push(boost::iostreams::gzip_decompressor());
	  inData[i].push(ifData[i]);
	}
	uint32_t allEOF = 0;
	std::vector<bool> eof(c.nTmpFile, false);
	std::vector<std::pair<unsigned, unsigned> > lastPair(c.nTmpFile, std::make_pair(0, 0));
	std::vector<bool> lastPairValid(c.nTmpFile, false);
	while(allEOF < c.nTmpFile) {
	  // Next sorted record
	  int32_t bestIdx = -1;
	  for(uint32_t i = 0; i < c.nTmpFile; ++i) {
	    if (!eof[i]) {
	      // Read next hash pair if necessary
	      if (!lastPairValid[i]) {
		// Read next line
		std::istream incoming(&inData[i]);
		std::string line;
		if (getline(incoming, line)) {
		  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
		  boost::char_separator<char> sep("\t");
		  Tokenizer tokens(line, sep);
		  Tokenizer::iterator tokIter = tokens.begin();
		  unsigned h1 = boost::lexical_cast<unsigned>(*tokIter++);
		  unsigned h2 = boost::lexical_cast<unsigned>(*tokIter++);
		  lastPair[i].first = h1;
		  lastPair[i].second = h2;
		  lastPairValid[i] = true;
		} else {
		  eof[i] = true;
		  ++allEOF;
		}
	      }
	      // Find next sorted hash pair
	      if (lastPairValid[i]) {
		if (bestIdx == -1) bestIdx = i;
		else {
		  if ((lastPair[i].first < lastPair[bestIdx].first) || ((lastPair[i].first == lastPair[bestIdx].first) && (lastPair[i].second < lastPair[bestIdx].second))) bestIdx = i;
		}
	      }
	    }
	  }
	  if (bestIdx != -1) {
	    dataOut << lastPair[bestIdx].first << '\t' << lastPair[bestIdx].second << std::endl;
	    lastPairValid[bestIdx] = false;
	  } else {
	    if (allEOF < c.nTmpFile) {
	      std::cerr << "Error: Unknown next sorted hash!" << std::endl;
	      exit(-1);
	    }
	  }
	}
	for(uint32_t i = 0; i < c.nTmpFile; ++i) {
	  inData[i].pop();
	  inData[i].pop();

	  // Delete sorted chunk
	  std::string read1sort = c.fq1 + "." + boost::lexical_cast<std::string>(i) + ".sort.gz";
	  boost::filesystem::remove(read1sort);
	}

	// Close sorted output file
	dataOut.pop();
	dataOut.pop();

	// Done
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Merged " << c.nTmpFile << " chunks." << std::endl;	
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
