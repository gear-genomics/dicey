#ifndef DESIGN_H
#define DESIGN_H

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

#include <sdsl/suffix_arrays.hpp>
#include <htslib/faidx.h>

#include "neighbors.h"
#include "align.h"
#include "needle.h"
#include "thal.h"
#include "util.h"
#include "gtf.h"
#include "gff3.h"

using namespace sdsl;

namespace dicey
{

  struct DesignConfig {
    bool indel;
    bool pruneprimer;
    bool hasOutfile;
    
    double cutTemp;
    uint32_t maxProdSize;
    double cutofPen;
    double penDiff;
    double penMis;
    double penLen;
    uint32_t kmer;
    uint32_t distance;
    uint32_t maxNeighborhood;
    uint32_t maxPruneCount;
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
    std::string idname;
    std::string feature;
    std::vector<std::string> chrname;
    std::map<std::string, int32_t> nchr;
    boost::filesystem::path gtfFile;
    boost::filesystem::path primer3Config;
    boost::filesystem::path outfile;
    boost::filesystem::path infile;
    boost::filesystem::path genome;
  };

  int design(int argc, char** argv) {
    DesignConfig c;
    
    // CMD Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("gtf,t", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3 file")
      ("config,i", boost::program_options::value<boost::filesystem::path>(&c.primer3Config)->default_value("./src/primer3_config/"), "primer3 config directory")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output file")
      ;

    boost::program_options::options_description appr("Approximate Search Options");
    appr.add_options()
      ("kmer,k", boost::program_options::value<uint32_t>(&c.kmer)->default_value(15), "k-mer size")
      ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(10000), "max. number of matches per k-mer")
      ("maxNeighborhood,x", boost::program_options::value<uint32_t>(&c.maxNeighborhood)->default_value(10000), "max. neighborhood size")
      ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
      ("pruneprimer,q", boost::program_options::value<uint32_t>(&c.maxPruneCount), "prune primer threshold")
      ("hamming,n", "use hamming neighborhood instead of edit distance")
    ;
    
    boost::program_options::options_description score("Parameters for Scoring and Penalty Calculation");
    score.add_options()
      ("cutTemp,c", boost::program_options::value<double>(&c.cutTemp)->default_value(45.0), "min. primer melting temperature")
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

    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Error: Genome does not exist!" << std::endl;
      return 1;
    }
    
    // Initialize thal arguments
    if ((!boost::filesystem::exists(c.primer3Config)) || (!boost::filesystem::is_directory(c.primer3Config))) {
      std::cerr << "Error: Cannot find primer3 config directory!" << std::endl;
      return 1;
    } else {
      c.primer3Config.remove_trailing_separator();
      c.primer3Config += boost::filesystem::path::preferred_separator;
      boost::filesystem::path filePath = c.primer3Config;
      filePath += "tetraloop.dh";
      if (!boost::filesystem::exists(filePath)) {
	std::cerr << "Error: Config directory path appears to be incorrect!" << std::endl;
	return 1;
      }
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
    // Fix provided Temperature
    a.temp += primer3thal::ABSOLUTE_ZERO;
            
    // Set prefix and suffix based on edit distance
    c.pre_context = 0;
    c.post_context = 0;
    if (c.indel) {
      c.pre_context += c.distance;
      c.post_context += c.distance;
    }

    // Cmd switches
    if (!vm.count("hamming")) c.indel = true;
    else c.indel = false;
    if (vm.count("outfile")) c.hasOutfile = true;
    else c.hasOutfile = false;
    if (vm.count("pruneprimer")) c.pruneprimer = true;
    else c.pruneprimer = false;
    
    // Fill genome map
    if (c.nchr.empty()) {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      c.chrname.resize(faidx_nseq(fai));
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string chrName = faidx_iseq(fai, refIndex);
	//std::cerr << chrName << ',' << refIndex << std::endl;
	c.nchr.insert(std::make_pair(chrName, refIndex));
	c.chrname[refIndex] = chrName;
      }
      fai_destroy(fai);
    }

    // Filter exon IDs
    c.idname = "gene_id";
    c.feature = "exon";
    
    // Parse exons
    typedef std::vector<IntervalLabel> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(c.nchr.size(), TChromosomeRegions());
    typedef std::vector<std::string> TGeneIds;
    TGeneIds geneIds;
    typedef std::vector<bool> TProteinCoding;
    TProteinCoding pCoding;
    parseGTF(c, gRegions, geneIds, pCoding);
    

    // Parse chromosomes
    faidx_t* fai = fai_load(c.genome.string().c_str());	
    for(uint32_t refIndex = 0; refIndex < c.nchr.size(); ++refIndex) {
      for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i) {
	if ((geneIds[gRegions[refIndex][i].lid] != "ENSG00000164692") && (geneIds[gRegions[refIndex][i].lid] != "ENSG00000172270")) continue;
	std::cerr << c.chrname[refIndex] << ':' << gRegions[refIndex][i].start << '-' << gRegions[refIndex][i].end << '\t' << gRegions[refIndex][i].strand << '\t' << gRegions[refIndex][i].lid << '\t' << geneIds[gRegions[refIndex][i].lid] << '\t' << pCoding[gRegions[refIndex][i].lid] << std::endl;

	// Parameters
	uint32_t armlen = 20;
	uint32_t targetlen = 2 * armlen;
	double armTMMax = 60;
	double armTMDiff = 2;
	double probeTMMin = 65;
	double probeTMMax = 75;
	double minGC = 0.4;
	double maxGC = 0.6;
	
	int32_t seqlen;
	char* seq = faidx_fetch_seq(fai, c.chrname[refIndex].c_str(), gRegions[refIndex][i].start, gRegions[refIndex][i].end, &seqlen);
	std::string exonseq = boost::to_upper_copy(std::string(seq));
	if (gRegions[refIndex][i].strand == '-') revcomplement(exonseq);
	uint32_t exonlen = exonseq.size();
	if (exonlen >= targetlen) {
	  std::string rexonseq(exonseq);
	  revcomplement(rexonseq);
	  for(uint32_t k = 0; k < (exonlen - targetlen + 1); ++k) {
	    // Arm1
	    std::string arm1 = exonseq.substr(k, armlen);
	    double arm1GC = gccontent(arm1);
	    if ((arm1GC < minGC) || (arm1GC > maxGC)) continue;
	    std::string rarm1 = rexonseq.substr(exonlen - armlen - k, armlen);
	    primer3thal::oligo1 = (unsigned char*) arm1.c_str();
	    primer3thal::oligo2 = (unsigned char*) rarm1.c_str();
	    primer3thal::thal_results oiarm1;
	    bool thalsuccess1 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oiarm1);
	    if ((!thalsuccess1) || (oiarm1.temp == primer3thal::THAL_ERROR_SCORE)) {
	      std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
	      return 1;
	    }
	    double arm1TM = oiarm1.temp;
	    if (arm1TM > armTMMax) continue;

	    // Arm2
	    std::string arm2 = exonseq.substr(k + armlen, armlen);
	    double arm2GC = gccontent(arm2);
	    if ((arm2GC < minGC) || (arm2GC > maxGC)) continue;
	    std::string rarm2 = rexonseq.substr(exonlen - armlen - (k + armlen), armlen);
	    primer3thal::oligo1 = (unsigned char*) arm2.c_str();
	    primer3thal::oligo2 = (unsigned char*) rarm2.c_str();
	    primer3thal::thal_results oiarm2;
	    bool thalsuccess2 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oiarm2);
	    if ((!thalsuccess2) || (oiarm2.temp == primer3thal::THAL_ERROR_SCORE)) {
	      std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
	      return 1;
	    }
	    double arm2TM = oiarm2.temp;
	    if ((arm2TM > armTMMax) || (std::abs(arm1TM - arm2TM) > armTMDiff)) continue;

	    // Probe
	    std::string probe = exonseq.substr(k, targetlen);
	    double probeGC = gccontent(probe);
	    if ((probeGC < minGC) || (probeGC > maxGC)) continue;
	    std::string rprobe = rexonseq.substr(exonlen - targetlen - k, targetlen);
	    primer3thal::oligo1 = (unsigned char*) probe.c_str();
	    primer3thal::oligo2 = (unsigned char*) rprobe.c_str();
	    primer3thal::thal_results oiprobe;
	    bool thalsuccess3 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oiprobe);
	    if ((!thalsuccess3) || (oiarm2.temp == primer3thal::THAL_ERROR_SCORE)) {
	      std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
	      return 1;
	    }
	    double probeTM = oiprobe.temp;
	    if ((probeTM < probeTMMin) || (probeTM > probeTMMax)) continue;

	    // Output
	    std::cerr << arm1 << '-' << arm2 << ",TM:" << arm1TM << ',' << arm2TM << ',' << probeTM << ",GC:" << arm1GC << ',' << arm2GC << ',' << probeGC << std::endl;
	  }
	}
	free(seq);
      }
    }
    fai_destroy(fai);
    exit(-1);
    


    // DNA Hits
    typedef std::vector<PrimerBind> TPrimerBinds;
    TPrimerBinds allp;
    typedef std::vector<PcrProduct> TPcrProducts;
    TPcrProducts pcrColl;
    std::vector<uint32_t> seqlen;
    std::vector<std::string> seqname;
    std::vector<std::string> pName;
    std::vector<std::string> pSeq;

    // Reference index
    csa_wt<> fm_index;  
    boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
    std::string index_file = op.string() + ".fm9";
    if (!load_from_checked_file(fm_index, index_file)) {
      std::cerr << "Error: FM-Index cannot be loaded!" << std::endl;
      return 1;
    }
    
    // Define alphabet
    typedef std::set<char> TAlphabet;
    char tmp[] = {'A', 'C', 'G', 'T'};
    TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
    
    // Query FM-Index
    typedef std::vector<TPrimerBinds> TChrPrimerBinds;
    TChrPrimerBinds forBind(c.chrname.size(), TPrimerBinds());
    TChrPrimerBinds revBind(c.chrname.size(), TPrimerBinds());
    for(uint32_t primerId = 0; primerId < pSeq.size(); ++primerId) {
      // Thermodynamic calculation
      std::string forQuery = pSeq[primerId];
      std::string revQuery = pSeq[primerId];
      reverseComplement(revQuery);
      primer3thal::oligo1 = (unsigned char*) forQuery.c_str();
      primer3thal::oligo2 = (unsigned char*) revQuery.c_str();
      primer3thal::thal_results oi;
      bool thalsuccess1 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oi);
      if ((!thalsuccess1) || (oi.temp == primer3thal::THAL_ERROR_SCORE)) {
	std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
	return 1;
      }
      double matchTemp = oi.temp;

      // Enumerate neighbors
      typedef std::set<std::string> TStringSet;
      typedef std::vector<TStringSet> TFwdRevSearchSets;
      TFwdRevSearchSets fwrv(2, TStringSet());
      std::string sequence = pSeq[primerId];
      uint32_t koffset = sequence.size() - c.kmer;
      sequence = sequence.substr(sequence.size() - c.kmer);
      neighbors(sequence, alphabet, c.distance, c.indel, c.maxNeighborhood, fwrv[0]);
      std::string revSequence = sequence;
      reverseComplement(revSequence);
      neighbors(revSequence, alphabet, c.distance, c.indel, c.maxNeighborhood, fwrv[1]);
      if ((fwrv[0].size() >= c.maxNeighborhood) || (fwrv[1].size() >= c.maxNeighborhood)) {
	std::cerr << "Warning: Neighborhood size exceeds " + boost::lexical_cast<std::string>(c.maxNeighborhood) + " candidates. Only first " + boost::lexical_cast<std::string>(c.maxNeighborhood) + " neighbors are searched, results are likely incomplete!" << std::endl;
      }

      // Serach
      uint32_t hits = 0;
      for(uint32_t fwrvidx = 0; fwrvidx < fwrv.size(); ++fwrvidx) {
	typedef std::pair<uint32_t, uint32_t> TRefPosPair;
	typedef std::set<TRefPosPair> TUniquePrimerHits;
	TUniquePrimerHits uphit;
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
	      if (fwrvidx) post_extract += koffset;
	      else pre_extract += koffset;
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
	      
	      // Thermodynamic calculation
	      std::string genomicseq = pre + s.substr(0, m) + post;
	      if (pre.size() < chrpos) chrpos -= pre.size();	      
	      std::string primer = revQuery;
	      std::string searchSeq = sequence;
	      if (fwrvidx) {
		primer = forQuery;
		searchSeq = revSequence;
	      }
	      primer3thal::oligo1 = (unsigned char*) primer.c_str();
	      primer3thal::oligo2 = (unsigned char*) genomicseq.c_str();
	      primer3thal::thal_results o;
	      bool thalsuccess = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &o);
	      if ((!thalsuccess) || (o.temp == primer3thal::THAL_ERROR_SCORE)) {
		std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
		return 1;
	      }

	      // Suitable match?
	      if (o.temp > c.cutTemp) {
		// Based on search sequence we can define the unique starting position to remove duplicates
		uint32_t alignpos = chrpos;
		DnaScore<int32_t> sc(0, -1, -1, -1);
		typedef boost::multi_array<char, 2> TAlign;
		AlignConfig<false, true> global;
		TAlign align;
		needle(genomicseq, searchSeq, align, global, sc);
		// Determine alignpos
		bool leadGap = true;
		for(uint32_t j = 0; (j < (align.shape()[1] - _trailGap(align))); ++j) {
		  if (align[1][j] != '-') leadGap = false;
		  if (leadGap) ++alignpos;
		}
		if (uphit.find(std::make_pair(refIndex, alignpos)) == uphit.end()) {
		  // New hit
		  uphit.insert(std::make_pair(refIndex, alignpos));

		  // Genomic subsequence
		  if (fwrvidx) {
		    chrpos = alignpos;
		    genomicseq = genomicseq.substr(alignpos - chrpos + 1, primer.size());
		  } else {
		    uint32_t alignshift = alignpos - chrpos;
		    chrpos = alignpos - koffset;
		    if (alignshift >= koffset) {
		      alignshift -= koffset;
		      genomicseq = genomicseq.substr(alignshift,  primer.size());
		    }
		  }

		  // New primer record
		  PrimerBind prim;
		  prim.refIndex = refIndex;
		  prim.temp = o.temp;
		  prim.perfTemp = matchTemp;
		  prim.primerId = primerId;
		  prim.genome = genomicseq;
		  if (fwrvidx) {
		    prim.onFor = false;
		    prim.pos = chrpos;
		    revBind[refIndex].push_back(prim);
		  } else {
		    prim.onFor = true;
		    prim.pos = chrpos;
		    forBind[refIndex].push_back(prim);
		  }
		}
	      }
	      ++hits;
	    }
	  }
	}
      }
      if (hits >= c.max_locations) {
	std::cerr << "Warning: More than " + boost::lexical_cast<std::string>(c.max_locations) + " matches found. Only first " + boost::lexical_cast<std::string>(c.max_locations) + " matches are reported, results are likely incomplete!" << std::endl;
      }
    }

    // Collect all primers
    for(uint32_t refIndex = 0; refIndex < c.chrname.size(); ++refIndex) {
      allp.insert( allp.end(), forBind[refIndex].begin(), forBind[refIndex].end() );
      allp.insert( allp.end(), revBind[refIndex].begin(), revBind[refIndex].end() );
    }
    
    // Sort by temperature
    std::sort(allp.begin(), allp.end(), SortPrimer<PrimerBind>());
    
    // Search PCR amplicons
    if (c.pruneprimer) {
      //jsonPrimerOut(c, seqname, allp, pcrColl, pName, pSeq, msg);
    } else {
      for(uint32_t refIndex = 0; refIndex < c.chrname.size(); ++refIndex) {
	for(TPrimerBinds::iterator fw = forBind[refIndex].begin(); fw != forBind[refIndex].end(); ++fw) {
	  for(TPrimerBinds::iterator rv = revBind[refIndex].begin(); rv != revBind[refIndex].end(); ++rv) {
	    if ((rv->pos > fw->pos) && (rv->pos - fw->pos < c.maxProdSize)) {
	      PcrProduct pcrProd;
	      pcrProd.refIndex = refIndex;
	      pcrProd.forPos = fw->pos;
	      pcrProd.forTemp = fw->temp;
	      pcrProd.forId = fw->primerId;
	      pcrProd.revPos = rv->pos;
	      pcrProd.revTemp = rv->temp;
	      pcrProd.revId = rv->primerId;
	      pcrProd.leng = (rv->pos + pSeq[pcrProd.revId].size()) - fw->pos;
	      
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
      //jsonPrimerOut(c, seqname, allp, pcrColl, pName, pSeq, msg);
    }
    
    // Clean-up
    primer3thal::destroy_thal_structures();
    
    return 0;
  }

}

#endif
