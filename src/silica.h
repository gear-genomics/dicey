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

#ifndef SILICA_H
#define SILICA_H

#include <iostream>
#include <fstream>
#include <vector>

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
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

#include <sdsl/suffix_arrays.hpp>
#include <htslib/faidx.h>

#include "neighbors.h"
#include "align.h"
#include "needle.h"
#include "thal.h"
#include "json.h"
#include "util.h"

using namespace sdsl;

namespace dicey
{

  struct SilicaConfig {
    bool indel;
    bool pruneprimer;
    double cutTemp;
    uint32_t maxProdSize;
    double cutofPen;
    double penDiff;
    double penMis;
    double penLen;
    uint32_t kmer;
    uint32_t distance;
    uint32_t maxNeighborhood;

    // Primer3
    double temp;
    double mv;
    double dv;
    double dna_conc;
    double dntp;
    std::size_t pre_context;
    std::size_t post_context;
    std::size_t max_locations;
    std::string format;
    boost::filesystem::path primer3Config;
    boost::filesystem::path outfile;
    boost::filesystem::path primfile;
    boost::filesystem::path infile;
    boost::filesystem::path genome;
  };
  
  struct PrimerBind {
    uint32_t refIndex;
    uint32_t pos;
    uint32_t primerId;
    bool onFor;
    double temp;
    double perfTemp;
    std::string genSeq;
  };
  
  template<typename TRecord>
    struct SortPrimer : public std::binary_function<TRecord, TRecord, bool>
    {
      inline bool operator()(TRecord const& a, TRecord const& b) const {
	return (a.temp > b.temp);
      }
    };

  struct PcrProduct {
    uint32_t refIndex;
    uint32_t leng;
    uint32_t forPos;
    uint32_t revPos;
    uint32_t forId;
    uint32_t revId;
    double forTemp;
    double revTemp;
    double penalty;
  };
  
  template<typename TRecord>
    struct SortProducts : public std::binary_function<TRecord, TRecord, bool>
    {
      inline bool operator()(TRecord const& a, TRecord const& b) const {
	return (a.penalty < b.penalty);
      }
    };
  
  template<typename TPrimerBinds>
    inline void
    addUnique(TPrimerBinds& coll, PrimerBind& prim, uint32_t const distance) {
    uint32_t idx = 0;
    for(typename TPrimerBinds::iterator it = coll.begin(); it != coll.end(); ++it, ++idx) {
      if ((prim.primerId == it->primerId) && (prim.pos + 2*distance >= it->pos) && (prim.pos <= it->pos + 2*distance)) {
	// Duplicate, find best temp
	if (prim.temp > it->temp) coll[idx] = prim;
	return;
      }
    }
    // No other primer found
    coll.push_back(prim);
  }
  
  
  int silica(int argc, char** argv) {
    SilicaConfig c;
    
    // CMD Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("config,i", boost::program_options::value<boost::filesystem::path>(&c.primer3Config)->default_value("./src/primer3_config/"), "primer3 config directory")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("amplicons.txt"), "amplicon output file")
      ("primer,p", boost::program_options::value<boost::filesystem::path>(&c.primfile)->default_value("primers.txt"), "primer locations file")
      ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("txt"), "output format (json, txt, csv or jsoncsv)")
      ;

    boost::program_options::options_description appr("Approximate Search Options");
    appr.add_options()
      ("kmer,k", boost::program_options::value<uint32_t>(&c.kmer)->default_value(15), "k-mer size")
      ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(10000), "max. number of matches per k-mer")
      ("maxNeighborhood,x", boost::program_options::value<uint32_t>(&c.maxNeighborhood)->default_value(10000), "max. neighborhood size")
      ("pruneprimer,q", "prune primers with more than maxmatches") 
      ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
      ("hamming,n", "use hamming neighborhood instead of edit distance")
    ;
    
    boost::program_options::options_description score("Parameters for Scoring and Penalty Calculation");
    score.add_options()
      ("cutTemp,c", boost::program_options::value<double>(&c.cutTemp)->default_value(40.0), "min. primer melting temperature")
      ("maxProdSize,l", boost::program_options::value<uint32_t>(&c.maxProdSize)->default_value(15000), "max. PCR Product size")
      ("cutoffPenalty", boost::program_options::value<double>(&c.cutofPen)->default_value(-1.0), "max. penalty for products (-1 = keep all)")
      ("penaltyTmDiff", boost::program_options::value<double>(&c.penDiff)->default_value(0.6), "multiplication factor for deviation of primer Tm penalty")
      ("penaltyTmMismatch", boost::program_options::value<double>(&c.penMis)->default_value(0.4), "multiplication factor for Tm pair difference penalty")
      ("penaltyLength", boost::program_options::value<double>(&c.penLen)->default_value(0.001), "multiplication factor for amplicon length penalty")
      ;

    boost::program_options::options_description tmcalc("Parameters for Tm Calculation");
    tmcalc.add_options()
      ("enttemp", boost::program_options::value<double>(&c.temp)->default_value(37.0), "temperature for entropie and entalpie calculation in Celsius")
      ("monovalent", boost::program_options::value<double>(&c.mv)->default_value(50.0), "concentration of monovalent ions in mMol")
      ("divalent", boost::program_options::value<double>(&c.dv)->default_value(1.5), "concentration of divalent ions in mMol")
      ("dna", boost::program_options::value<double>(&c.dna_conc)->default_value(50.0), "concentration of annealing(!) Oligos in nMol")
      ("dntp", boost::program_options::value<double>(&c.dntp)->default_value(0.6), "the sum  of all dNTPs in mMol")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "seq.fasta")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(appr).add(score).add(tmcalc).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(appr).add(score).add(tmcalc);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] -g <ref.fa.gz> sequences.fasta" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    
    // Cmd switches
    if (!vm.count("hamming")) c.indel = true;
    else c.indel = false;
    if (vm.count("pruneprimer")) c.pruneprimer = true;
    else c.pruneprimer = false;

    // Initialize thal arguments
    if (!boost::filesystem::exists(c.primer3Config)) {
      std::cerr << "Cannot find primer3 config directory: " << c.primer3Config.string() << std::endl;
      return 1;
    }
    primer3thal::thal_args a;
    primer3thal::set_thal_default_args(&a);
    a.temponly=1;
    a.type = primer3thal::thal_end1;
    primer3thal::get_thermodynamic_values(c.primer3Config.string().c_str());
    a.temp = c.temp;
    a.mv = c.mv;
    a.dv = c.dv;
    a.dna_conc = c.dna_conc;
    a.dntp = c.dntp;
    
    // Set prefix and suffix based on edit distance
    c.pre_context = 2;
    c.post_context = 2;
    if (c.indel) {
      c.pre_context += c.distance;
      c.post_context += c.distance;
    }
    
    // Fix provided Temperature
    a.temp += primer3thal::ABSOLUTE_ZERO;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "dicey ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parse chromosomes" << std::endl;
    
    // Parse chromosome lengths
    std::vector<uint32_t> seqlen;
    std::vector<std::string> seqname;
    uint32_t nseq = getSeqLenName(c, seqlen, seqname);
    if (!nseq) {
      std::cerr << "Could not retrieve sequence lengths!" << std::endl;
      return 1;
    }

    // Open fasta index
    faidx_t* fai = fai_load(c.genome.string().c_str());
  
    // Reference index
    csa_wt<> fm_index;  
    
    // Load FM index
    boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
    std::string index_file = op.string() + ".fm9";
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FM-Index" << std::endl;
    if (!load_from_file(fm_index, index_file)) {
      std::cerr << "Index cannot be loaded!" << std::endl;
      return 1;
    }
    
    // Parse input fasta
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parse input FASTA" << std::endl;
    
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      std::cerr << "Input fasta file is missing: " << c.infile.string() << std::endl;
      return 1;
    }
    typedef std::vector<std::string> TPrimerName;
    typedef std::vector<std::string> TPrimerSeq;
    TPrimerName pName;
    TPrimerSeq pSeq;
    std::ifstream fafile(c.infile.string().c_str());
    if (fafile.good()) {
      std::string fan;
      std::string tmpfasta;
      std::string line;
      while(std::getline(fafile, line)) {
	if (!line.empty()) {
	  if (line[0] == '>') {
	    if ((!fan.empty()) && (!tmpfasta.empty()) && (tmpfasta.size() > c.kmer)) {
	      if (c.pruneprimer) {
		std::string qr = tmpfasta.substr(tmpfasta.size() - c.kmer);
		std::size_t occs = sdsl::count(fm_index, qr.begin(), qr.end());
		if (occs <= c.max_locations) {
		  reverseComplement(qr);
		  occs = sdsl::count(fm_index, qr.begin(), qr.end());
		  if (occs <= c.max_locations) {
		    pName.push_back(fan);
		    pSeq.push_back(tmpfasta);
		  }
		}
	      } else {
		pName.push_back(fan);
		pSeq.push_back(tmpfasta);
	      }
	      tmpfasta = "";
	    }
	    fan = line.substr(1);
	  } else {
	    tmpfasta += boost::to_upper_copy(line);
	  }
	}
      }
      if ((!fan.empty()) && (!tmpfasta.empty()) && (tmpfasta.size() > c.kmer)) {
	if (c.pruneprimer) {
	  std::string qr = tmpfasta.substr(tmpfasta.size() - c.kmer);
	  std::size_t occs = sdsl::count(fm_index, qr.begin(), qr.end());
	  if (occs <= c.max_locations) {
	    reverseComplement(qr);
	    occs = sdsl::count(fm_index, qr.begin(), qr.end());
	    if (occs <= c.max_locations) {
	      pName.push_back(fan);
	      pSeq.push_back(tmpfasta);
	    }
	  }
	} else {
	  pName.push_back(fan);
	  pSeq.push_back(tmpfasta);
	}
      }
    }
    
    // Define alphabet
    typedef std::set<char> TAlphabet;
    char tmp[] = {'A', 'C', 'G', 'T'};
    TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
    
    // Query FM-Index
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Query FM-Index" << std::endl;
    boost::progress_display show_progress( pSeq.size() );
    typedef std::vector<PrimerBind> TPrimerBinds;
    typedef std::vector<TPrimerBinds> TChrPrimerBinds;
    TChrPrimerBinds forBind(faidx_nseq(fai), TPrimerBinds());
    TChrPrimerBinds revBind(faidx_nseq(fai), TPrimerBinds());
    for(uint32_t primerId = 0; primerId < pSeq.size(); ++primerId) {
      std::string qr = pSeq[primerId];
      if (qr.size() < c.kmer) continue;
      int32_t koffset = qr.size() - c.kmer;
      qr = qr.substr(qr.size() - c.kmer);
      typedef std::set<std::string> TStringSet;
      TStringSet fwdset;
      neighbors(qr, alphabet, c.distance, c.indel, c.maxNeighborhood, fwdset);
      // Debug
      //for(TStringSet::iterator it = fwdset.begin(); it != fwdset.end(); ++it) std::cerr << *it << std::endl;
      TStringSet revset;
      reverseComplement(qr);
      neighbors(qr, alphabet, c.distance, c.indel, c.maxNeighborhood, revset);
      int32_t qhits = 0;
      for(int32_t fwdrev = 0; fwdrev < 2; ++fwdrev) {
	TStringSet::iterator its;
	TStringSet::iterator itsEnd;
	if (fwdrev == 0) {
	  its = fwdset.begin();
	  itsEnd = fwdset.end();
	} else {
	  its = revset.begin();
	  itsEnd = revset.end();
	}
	for(; its != itsEnd; ++its, ++qhits) {
	  std::string query(*its);
	  std::string forQuery = pSeq[primerId];
	  std::string revQuery = pSeq[primerId];
	  reverseComplement(revQuery);
	  primer3thal::oligo1 = (unsigned char*) forQuery.c_str();
	  primer3thal::oligo2 = (unsigned char*) revQuery.c_str();
	  primer3thal::thal_results oi;
	  bool thalsuccess1 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oi);
	  if ((!thalsuccess1) || (oi.temp == primer3thal::THAL_ERROR_SCORE)) {
	    std::cerr << "Error during thermodynamical calculation!" << std::endl;
	    return -1;
	  }
	  double matchTemp = oi.temp;
	  std::size_t m = query.size();
	  std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
	  if (occs > 0) {
	    auto locations = locate(fm_index, query.begin(), query.begin() + m);
	    std::sort(locations.begin(), locations.end());
	    std::size_t pre_extract = c.pre_context;
	    std::size_t post_extract = c.post_context;
	    if (fwdrev == 0) pre_extract += koffset;
	    else post_extract += koffset;
	    for(std::size_t i = 0; i < std::min(occs, c.max_locations); ++i) {
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
	      
	      // thermodynamical alignemnt
	      std::string genomicseq = pre + s.substr(0, m) + post;
	      std::string primer = pSeq[primerId];
	      if (fwdrev == 0) reverseComplement(primer);
	      primer3thal::oligo1 = (unsigned char*) primer.c_str();
	      primer3thal::oligo2 = (unsigned char*) genomicseq.c_str();
	      primer3thal::thal_results o;
	      bool thalsuccess = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &o);
	      if ((!thalsuccess) || (o.temp == primer3thal::THAL_ERROR_SCORE)) {
		std::cerr << "Error during thermodynamical calculation!" << std::endl;
		return -1;
	      }
	      
	      // Score suitable primers
	      if (o.temp > c.cutTemp) {
		PrimerBind prim;
		prim.refIndex = refIndex;
		prim.temp = o.temp;
		prim.perfTemp = matchTemp;
		prim.primerId = primerId;
		prim.genSeq = genomicseq;
		if (fwdrev == 0) {
		  prim.onFor = true;
		  prim.pos = chrpos - koffset;
		  if (c.indel) {
		    prim.pos -= c.distance;
		    addUnique(forBind[refIndex], prim, c.distance);
		  } else forBind[refIndex].push_back(prim);
		} else {
		  prim.onFor = false;
		  prim.pos = chrpos + pSeq[primerId].size();
		  if (c.indel) {
		    prim.pos += c.distance;
		    addUnique(revBind[refIndex], prim, c.distance);
		  } else revBind[refIndex].push_back(prim);
		}
	      }
	      
	      // Debug alignment output
	      //_debugAlignment(pSeq[primerId], genomicseq, fwdrev);
	    }
	  }
	}
      }
      ++show_progress;
    }

    // Collect all primers
    TPrimerBinds allp;  
    for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
      allp.insert( allp.end(), forBind[refIndex].begin(), forBind[refIndex].end() );
      allp.insert( allp.end(), revBind[refIndex].begin(), revBind[refIndex].end() );
    }
    
    // Sort by temperature
    std::sort(allp.begin(), allp.end(), SortPrimer<PrimerBind>());
    
    // Output primers
    if (c.format == "json") primerJsonOut(c.primfile.string(), fai, allp, pName, pSeq);
    else if (c.format == "csv") primerCsvOut(c.primfile.string(), fai, allp, pName, pSeq);
    else if (c.format == "jsoncsv") {
      primerJsonOut(c.primfile.string() + ".json", fai, allp, pName, pSeq);
      primerCsvOut(c.primfile.string() + ".csv", fai, allp, pName, pSeq);
    }
    else primerTxtOut(c.primfile.string(), fai, allp, pName, pSeq);
    
    if (!c.pruneprimer) {
      // Find PCR products
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Find PCR Products" << std::endl;
      boost::progress_display sp( faidx_nseq(fai) );
      typedef std::vector<PcrProduct> TPcrProducts;
      TPcrProducts pcrColl;
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	++sp;
	for(TPrimerBinds::iterator fw = forBind[refIndex].begin(); fw != forBind[refIndex].end(); ++fw) {
	  for(TPrimerBinds::iterator rv = revBind[refIndex].begin(); rv != revBind[refIndex].end(); ++rv) {
	    if ((rv->pos > fw->pos) && (rv->pos - fw->pos < c.maxProdSize)) {
	      PcrProduct pcrProd;
	      pcrProd.leng = rv->pos - fw->pos;
	      pcrProd.refIndex = refIndex;
	      pcrProd.forPos = fw->pos;
	      pcrProd.forTemp = fw->temp;
	      pcrProd.forId = fw->primerId;
	      pcrProd.revPos = rv->pos;
	      pcrProd.revTemp = rv->temp;
	      pcrProd.revId = rv->primerId;
	      // Calculate Penalty
	      double pen = (fw->perfTemp - fw->temp) * c.penDiff;
	      if (pen < 0) pen = 0;
	      double bpen = (rv->perfTemp - rv->temp) * c.penDiff;
	      if (bpen > 0) pen += bpen;
	      pen += std::abs(fw->temp - rv->temp) * c.penMis;
	      pen += pcrProd.leng * c.penLen;
	      pcrProd.penalty = pen;
	      if ((c.cutofPen < 0) || (pen < c.cutofPen)) pcrColl.push_back(pcrProd);
	    }
	  }
      }
      }
      
      // Sort by penalty
      std::sort(pcrColl.begin(), pcrColl.end(), SortProducts<PcrProduct>());
      
      // Output amplicons
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output Amplicons" << std::endl;
      if (c.format == "json") ampliconJsonOut(c.outfile.string(), fai, pcrColl, pName, pSeq);
      else if (c.format == "csv") ampliconCsvOut(c.outfile.string(), fai, pcrColl, pName, pSeq);
      else if (c.format == "jsoncsv") {
	ampliconJsonOut(c.outfile.string() + ".json", fai, pcrColl, pName, pSeq);
      ampliconCsvOut(c.outfile.string() + ".csv", fai, pcrColl, pName, pSeq);
      }
      else ampliconTxtOut(c.outfile.string(), fai, pcrColl, pName, pSeq);
    }
    
    // Clean-up
    primer3thal::destroy_thal_structures();
    if (fai != NULL) fai_destroy(fai);
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

    return 0;
  }

}

#endif
