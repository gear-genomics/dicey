#ifndef PADLOCK_H
#define PADLOCK_H

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

  struct PadlockConfig {
    bool indel;
    bool nonprotein;
    bool computeAll;
    uint32_t distance;
    double temp;
    double mv;
    double dv;
    double dna_conc;
    double dntp;
    std::string idname;
    std::string feature;
    std::vector<std::string> chrname;
    std::map<std::string, int32_t> nchr;
    boost::filesystem::path gtfFile;
    boost::filesystem::path barcodes;
    boost::filesystem::path primer3Config;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path infile;
  };

  int padlock(int argc, char** argv) {
    PadlockConfig c;
    
    // CMD Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
      ("gtf,t", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3 file")
      ("barcodes,b", boost::program_options::value<boost::filesystem::path>(&c.barcodes), "FASTA barcode file")
      ("config,i", boost::program_options::value<boost::filesystem::path>(&c.primer3Config)->default_value("./src/primer3_config/"), "primer3 config directory")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.tsv"), "output file")
      ("hamming,n", "use hamming neighborhood instead of edit distance")
      ("nonprotein,p", "include non-protein coding genes")
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
      ("infile", boost::program_options::value<boost::filesystem::path>(&c.infile), "gene.lst")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("infile", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(tmcalc).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(tmcalc);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("infile")) || (!vm.count("genome")) || (!vm.count("barcodes")) || (!vm.count("gtf"))) {
      std::cout << "Usage: dicey " << argv[0] << " [OPTIONS] -g <ref.fa.gz> -t <ref.gtf.gz> -b <barcodes.fa.gz> <gene.list.file>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Cmd switches
    if (vm.count("hamming")) c.indel = false;
    else c.indel = true;
    if (vm.count("nonprotein")) c.nonprotein = true;
    else c.nonprotein = false;
    
    // Parameters
    uint32_t maxNeighborHits = 0;
    if (c.indel) maxNeighborHits = 2 * c.distance; // first base deletion and last base deletion give a hit 
    double armTMMax = 60;
    double armTMDiff = 2;
    double probeTMMin = 65;
    double probeTMMax = 75;
    double minGC = 0.4;
    double maxGC = 0.6;
    std::string anchorseq = "TGCGTCTATTTAGTGGAGCC";
    std::string spacerleft = "TCCTC";
    std::string spacerright = "TCTTT";
    
    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Error: Genome does not exist!" << std::endl;
      return 1;
    }

    // Check gene list
    std::set<std::string> geneset;
    c.computeAll = false;
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      if (c.infile.string() == "all") c.computeAll = true;
      else {
	std::cerr << "Error: Gene list does not exist!" << std::endl;
	return 1;
      }
    } else {
      std::ifstream geneFile(c.infile.string().c_str(), std::ifstream::in);
      if (geneFile.is_open()) {
        while (geneFile.good()) {
          std::string gline;
          getline(geneFile, gline);
          typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
          boost::char_separator<char> sep(" \t,;");
          Tokenizer tokens(gline, sep);
          Tokenizer::iterator tokIter = tokens.begin();
          if (tokIter!=tokens.end()) {
            std::string geneName=*tokIter++;
	    geneset.insert(geneName);
	  }
	}
	geneFile.close();
      }
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

    // Load barcodes
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Load barcodes" << std::endl;
    typedef std::vector<std::string> TBarcodes;
    TBarcodes barcodes(geneIds.size(), "NNNNNNNNNNNNNNNNNNNN");
    uint32_t numBarcodes = 0;
    if (!numBarcodes) {
      faidx_t* fai = fai_load(c.barcodes.string().c_str());
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	int32_t seqlen;
	char* seq = faidx_fetch_seq(fai, faidx_iseq(fai, refIndex), 0, 20, &seqlen);
	if (seqlen != 20) {
	  std::cerr << "Invalid barcode length! Should be 20bp." << std::endl;
	  return 1;
	}
	while ((numBarcodes < barcodes.size()) && (!c.nonprotein) && (!pCoding[numBarcodes])) ++numBarcodes;
	while ((numBarcodes < barcodes.size()) && (!c.computeAll) && (geneset.find(geneIds[numBarcodes]) == geneset.end())) ++numBarcodes;
	if (numBarcodes < barcodes.size()) barcodes[numBarcodes] = boost::to_upper_copy(std::string(seq));
	++numBarcodes;
	free(seq);
      }
      fai_destroy(fai);
    }
    // Debug
    //for(uint32_t i = 0; i < pCoding.size(); ++i) std::cerr << geneIds[i] << ',' << pCoding[i] << ',' << barcodes[i] << std::endl;
    
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

    // Outfile
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Compute padlocks" << std::endl;
    std::ofstream ofile(c.outfile.string().c_str());
    ofile << "Gene\tStrand\tExonCoordinates\tProbeSeq\tSpacerLeft\tAnchorSeq\tBarcodeSeq\tSpacerRight\tPadlockSeq\tArm1TM\tArm2TM\tBarcodeTM\tProbeTM\tArm1GC\tArm2GC\tBarcodeGC\tProbeGC" << std::endl;
    // Parse chromosomes
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(uint32_t armlen = 20; armlen < 21; ++armlen) {
      uint32_t targetlen = 2 * armlen;
      for(uint32_t refIndex = 0; refIndex < c.nchr.size(); ++refIndex) {
	for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i) {
	  //std::cerr << c.chrname[refIndex] << '\t' << gRegions[refIndex][i].start << '\t' << gRegions[refIndex][i].end << '\t' << gRegions[refIndex][i].strand << '\t' << geneIds[gRegions[refIndex][i].lid] << std::endl;
	  if ((c.computeAll) && (!c.nonprotein) && (!pCoding[gRegions[refIndex][i].lid])) continue;
	  if ((!c.computeAll) && (geneset.find(geneIds[gRegions[refIndex][i].lid]) == geneset.end())) continue;

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
	      
	      // Search neighbors of arm1
	      typedef std::set<std::string> TStringSet;
	      typedef std::vector<TStringSet> TFwdRevSearchSets;
	      TFwdRevSearchSets fwrv(2, TStringSet());
	      neighbors(arm1, alphabet, c.distance, c.indel, 10000, fwrv[0]);  // X-bp difference, at most 10000 neighbors
	      neighbors(rarm1, alphabet, c.distance, c.indel, 10000, fwrv[1]);
	      std::vector<uint32_t> hits(2, 0);
	      for(uint32_t fwrvidx = 0; fwrvidx < fwrv.size(); ++fwrvidx) {
		for(typename TStringSet::const_iterator it = fwrv[fwrvidx].begin(); ((it != fwrv[fwrvidx].end()) && (hits[0] + hits[1] <= maxNeighborHits)); ++it) {
		  std::string query = *it;
		  std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
		  hits[fwrvidx] += occs;
		}
	      }
	      if (hits[0] + hits[1] > maxNeighborHits) continue;
	      
	      // Search neighbors of arm2
	      fwrv.clear();
	      fwrv.resize(2, TStringSet());
	      neighbors(arm2, alphabet, c.distance, c.indel, 10000, fwrv[0]);  // 1-bp difference, at most 10000 neighbors
	      neighbors(rarm2, alphabet, c.distance, c.indel, 10000, fwrv[1]);
	      hits.clear();
	      hits.resize(2, 0);
	      for(uint32_t fwrvidx = 0; fwrvidx < fwrv.size(); ++fwrvidx) {
		for(typename TStringSet::const_iterator it = fwrv[fwrvidx].begin(); ((it != fwrv[fwrvidx].end()) && (hits[0] + hits[1] <= maxNeighborHits)); ++it) {
		  std::string query = *it;
		  std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
		  hits[fwrvidx] += occs;
		}
	      }
	      if (hits[0] + hits[1] > maxNeighborHits) continue;

	      // Final padlock
	      std::string padlock = rarm1 + spacerleft + anchorseq + barcodes[gRegions[refIndex][i].lid] + spacerright + rarm2;
	      double padlockGC = gccontent(padlock);
	      if ((padlockGC < minGC) || (padlockGC > maxGC)) continue;

	      // Barcode
	      std::string bartmp = barcodes[gRegions[refIndex][i].lid];
	      double barGC = gccontent(bartmp);
	      std::string rbartmp(bartmp);
	      revcomplement(rbartmp);
	      primer3thal::oligo1 = (unsigned char*) bartmp.c_str();
	      primer3thal::oligo2 = (unsigned char*) rbartmp.c_str();
	      primer3thal::thal_results oibartmp;
	      bool thalsuccess5 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oibartmp);
	      if ((!thalsuccess5) || (oibartmp.temp == primer3thal::THAL_ERROR_SCORE)) {
		std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
		return 1;
	      }
	      double barTM = oibartmp.temp;
	      
	      // Output
	      ofile << geneIds[gRegions[refIndex][i].lid] << '\t' << gRegions[refIndex][i].strand << '\t' << c.chrname[refIndex] << ':' << gRegions[refIndex][i].start << '-' << gRegions[refIndex][i].end << '\t';
	      ofile << arm1 << '-' << arm2 << '\t';
	      ofile << spacerleft << '\t' << anchorseq << '\t' << barcodes[gRegions[refIndex][i].lid] << '\t' << spacerright << '\t';
	      ofile << padlock << '\t';
	      ofile << arm1TM << '\t' << arm2TM << '\t' << barTM << '\t' << probeTM << '\t';
	      ofile << arm1GC << '\t' << arm2GC << '\t' << barGC << '\t' << probeGC << std::endl;
	    }
	  }
	  free(seq);
	}
      }
    }
    fai_destroy(fai);
    ofile.close();
    
    // Clean-up
    primer3thal::destroy_thal_structures();
    
    return 0;
  }

}

#endif
