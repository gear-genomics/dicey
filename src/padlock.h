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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <sdsl/suffix_arrays.hpp>
#include <htslib/faidx.h>

#include "neighbors.h"
#include "align.h"
#include "needle.h"
#include "thal.h"
#include "util.h"
#include "gtf.h"

using namespace sdsl;

namespace dicey
{

  struct PadlockConfig {
    bool json;
    bool indel;
    bool armMode;
    bool overlapping;
    bool computeAll;
    uint32_t distance;
    uint32_t armlen;
    double temp;
    double mv;
    double dv;
    double dna_conc;
    double dntp;
    std::string anchor;
    std::string spacerleft;
    std::string spacerright;    
    std::set<std::string> geneset;
    std::vector<std::string> chrname;
    std::map<std::string, int32_t> nchr;
    boost::filesystem::path gtfFile;
    boost::filesystem::path barcodes;
    boost::filesystem::path primer3Config;
    boost::filesystem::path outfile;
    boost::filesystem::path jsonfile;
    boost::filesystem::path genome;
    boost::filesystem::path infile;
  };


  inline void
  _outputGeneCodes() {
    // Generate all possible codes
    std::vector<std::string> geneCodes;
    std::vector<char> genealph;
    for(uint32_t i = 1; i < 5; ++i) genealph.push_back(boost::lexical_cast<char>(i));
    generateGeneCode(genealph, "", 6, geneCodes);

    // Pick at random
    boost::random::mt19937 gen;
    boost::random::uniform_int_distribution<> dist(0, geneCodes.size() - 1);
    boost::random::uniform_int_distribution<> deplete(0, 100);
    uint32_t maxCodes = 500;
    uint32_t minham = 2;
    std::vector<std::string> selected;
    uint32_t iterations = 0;
    while (selected.size() < maxCodes) {
      ++iterations;
      if (iterations % 10000000 == 0) std::cerr << "Iteration: " << iterations << std::endl;
      int32_t idx = dist(gen);
      std::set<char> diversity;
      int32_t onecount = 0;
      for(uint32_t k = 0; k < geneCodes[idx].size(); ++k) {
	diversity.insert(geneCodes[idx][k]);
	if (geneCodes[idx][k] == '1') ++onecount;
      }
      // Check diversity
      if (diversity.size() < 3) continue;
      // Deplete 1
      if (deplete(gen) > 100 - onecount * 49) continue;
      // Check hamming
      bool validham = true;
      for(uint32_t i = 0; i < selected.size(); ++i) {
	uint32_t ham = 0;
	for(uint32_t k = 0; k < geneCodes[idx].size(); ++k) {
	  if (geneCodes[idx][k] != selected[i][k]) ++ham;
	}
	if (ham < minham) {
	  validham = false;
	  break;
	}
      }
      if (validham) selected.push_back(geneCodes[idx]);
    }

    for(uint32_t i = 0; i < selected.size(); ++i) {
      std::cout << selected[i] << std::endl;
    }
  }
      

  template<typename TConfig>
  inline int32_t
  runPadlock(TConfig& c) {
#ifdef PROFILE
    ProfilerStart("dicey.prof");
#endif

    // Parameters
    uint32_t maxNeighborHits = 0;
    if (c.indel) maxNeighborHits = 2 * c.distance; // first base deletion and last base deletion give a hit 
    double armTMDiff = 2;
    double minGC = 0.4;
    double maxGC = 0.6;
    
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
            
    // Parse exons
    typedef std::vector<IntervalLabel> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(c.nchr.size(), TChromosomeRegions());
    typedef std::vector<GeneInfo> TGeneInfo;
    TGeneInfo geneInfo;
    parseGTF(c, gRegions, geneInfo);
    
    // Load barcodes
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Load barcodes" << std::endl;
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
	if (numBarcodes < geneInfo.size()) {
	  geneInfo[numBarcodes].barcode = boost::to_upper_copy(std::string(seq));
	  geneInfo[numBarcodes].code = std::string(faidx_iseq(fai, refIndex));
	}
	++numBarcodes;
	free(seq);
      }
      fai_destroy(fai);
    }
    
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

    // JSON output
    boost::iostreams::filtering_ostream rcfile;
    bool firstRec = true;
    if (c.json) {
      rcfile.push(boost::iostreams::gzip_compressor());
      rcfile.push(boost::iostreams::file_sink(c.jsonfile.c_str(), std::ios_base::out | std::ios_base::binary));
      rcfile << "{";
      rcfile << "\"errors\": [],";
      // Meta information
      rcfile << "\"meta\":";
      nlohmann::json meta;
      meta["version"] = diceyVersionNumber;
      meta["subcommand"] = "padlock";
      meta["armlength"] = c.armlen;
      meta["distance"] = c.distance;
      meta["distance"] = c.distance;
      meta["genome"] = c.genome.string();
      meta["infile"] = c.infile.string();
      meta["outfile"] = c.outfile.string();
      meta["barcodes"] = c.barcodes.string();
      meta["gtf"] = c.gtfFile.string();
      meta["jsonfile"] = c.jsonfile.string();
      meta["hamming"] = (!c.indel);
      rcfile << meta.dump() << ',';
      rcfile << "\"data\":{";
      rcfile << "\"columns\": [";
      rcfile << "\"Gene\", \"Symbol\", \"Code\", \"Position\", \"UCSC\", \"Strand\", \"ExonCoordinates\", \"ProbeSeq\", \"SpacerLeft\", \"AnchorSeq\", \"BarcodeSeq\", \"SpacerRight\", \"PadlockSeq\", \"Arm1TM\", \"Arm2TM\", \"BarcodeTM\", \"ProbeTM\", \"Arm1GC\", \"Arm2GC\", \"BarcodeGC\", \"ProbeGC\"";
      rcfile << "]," << std::endl;
      rcfile << "\"rows\": [" << std::endl;
    }
    
    // Outfile
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Compute padlocks" << std::endl;
    std::ofstream ofile(c.outfile.string().c_str());
    ofile << "Gene\tSymbol\tCode\tPosition\tUCSC\tStrand\tExonCoordinates\tProbeSeq\tSpacerLeft\tAnchorSeq\tBarcodeSeq\tSpacerRight\tPadlockSeq\tArm1TM\tArm2TM\tBarcodeTM\tProbeTM\tArm1GC\tArm2GC\tBarcodeGC\tProbeGC" << std::endl;
    // Parse chromosomes
    faidx_t* fai = fai_load(c.genome.string().c_str());
    uint32_t targetlen = 2 * c.armlen;
    for(uint32_t refIndex = 0; refIndex < c.nchr.size(); ++refIndex) {
      for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i) {
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
	    std::string arm1 = exonseq.substr(k, c.armlen);
	    double arm1GC = gccontent(arm1);
	    if ((arm1GC < minGC) || (arm1GC > maxGC)) continue;
	    std::string rarm1 = rexonseq.substr(exonlen - c.armlen - k, c.armlen);
	    primer3thal::oligo1 = (unsigned char*) arm1.c_str();
	    primer3thal::oligo2 = (unsigned char*) rarm1.c_str();
	    primer3thal::thal_results oiarm1;
	    bool thalsuccess1 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oiarm1);
	    if ((!thalsuccess1) || (oiarm1.temp == primer3thal::THAL_ERROR_SCORE)) {
	      std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
	      return 1;
	    }
	    double arm1TM = oiarm1.temp;
	    double armTMMax = 93 + arm1GC - 675 / c.armlen;  // ideal 81.5 instead of 93
	    if (arm1TM > armTMMax) continue;
	    
	    // Arm2
	    std::string arm2 = exonseq.substr(k + c.armlen, c.armlen);
	    double arm2GC = gccontent(arm2);
	    if ((arm2GC < minGC) || (arm2GC > maxGC)) continue;
	    std::string rarm2 = rexonseq.substr(exonlen - c.armlen - (k + c.armlen), c.armlen);
	    primer3thal::oligo1 = (unsigned char*) arm2.c_str();
	    primer3thal::oligo2 = (unsigned char*) rarm2.c_str();
	    primer3thal::thal_results oiarm2;
	    bool thalsuccess2 = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &oiarm2);
	    if ((!thalsuccess2) || (oiarm2.temp == primer3thal::THAL_ERROR_SCORE)) {
	      std::cerr << "Error: Thermodynamical calculation failed!" << std::endl;
	      return 1;
	    }
	    double arm2TM = oiarm2.temp;
	    armTMMax = 93 + arm2GC - 675 / c.armlen;
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
	    double probeTMMin = 81.5 + probeGC - 675 / (2 * c.armlen);
	    double probeTMMax = probeTMMin + 10;    
	    if ((probeTM < probeTMMin) || (probeTM > probeTMMax)) continue;

	    // Make sure the arms are unique
	    std::size_t ucount1 = sdsl::count(fm_index, arm1.begin(), arm1.end());
	    ucount1 += sdsl::count(fm_index, rarm1.begin(), rarm1.end());
	    if ((c.armMode) && (ucount1 > 1)) continue;
	    std::size_t ucount2 = sdsl::count(fm_index, arm2.begin(), arm2.end());
	    ucount2 += sdsl::count(fm_index, rarm2.begin(), rarm2.end());
	    if ((c.armMode) && (ucount2 > 1)) continue;

	    // Probe mode: At least one arm unique
	    if ((!c.armMode) && (ucount1 > 1) && (ucount2 > 1)) continue;
	    
	    // Search neighbors of arm1
	    typedef std::set<std::string> TStringSet;
	    typedef std::vector<TStringSet> TFwdRevSearchSets;
	    TFwdRevSearchSets fwrv(2, TStringSet());
	    neighbors(arm1, alphabet, c.distance, c.indel, 10000, fwrv[0]);  // X-bp difference, at most 10000 neighbors
	    neighbors(rarm1, alphabet, c.distance, c.indel, 10000, fwrv[1]);
	    
	    std::vector<uint32_t> hits(2, 0);
	    //std::cerr << gRegions[refIndex][i].strand << ',' << arm1 << ',' << rarm1 << std::endl;
	    for(uint32_t fwrvidx = 0; fwrvidx < fwrv.size(); ++fwrvidx) {
	      for(typename TStringSet::const_iterator it = fwrv[fwrvidx].begin(); ((it != fwrv[fwrvidx].end()) && (hits[0] + hits[1] <= maxNeighborHits)); ++it) {
		std::string query = *it;
		std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
		//if (occs > 0) std::cerr << *it << ',' << occs << std::endl;
		hits[fwrvidx] += occs;
	      }
	    }
	    if ((c.armMode) && (hits[0] + hits[1] > maxNeighborHits)) continue;
	    
	    // Search neighbors of arm2
	    fwrv.clear();
	    fwrv.resize(2, TStringSet());
	    neighbors(arm2, alphabet, c.distance, c.indel, 10000, fwrv[0]);  // 1-bp difference, at most 10000 neighbors
	    neighbors(rarm2, alphabet, c.distance, c.indel, 10000, fwrv[1]);
	    std::vector<uint32_t> hitsOther(2, 0);
	    for(uint32_t fwrvidx = 0; fwrvidx < fwrv.size(); ++fwrvidx) {
		for(typename TStringSet::const_iterator it = fwrv[fwrvidx].begin(); ((it != fwrv[fwrvidx].end()) && (hitsOther[0] + hitsOther[1] <= maxNeighborHits)); ++it) {
		  std::string query = *it;
		  std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
		  hitsOther[fwrvidx] += occs;
		}
	      }
	    if ((c.armMode) && (hitsOther[0] + hitsOther[1] > maxNeighborHits)) continue;

	    // Probe mode
	    if ((!c.armMode) && (hits[0] + hits[1] > maxNeighborHits) && (hitsOther[0] + hitsOther[1] > maxNeighborHits)) continue;
	    
	    // Final padlock
	    std::string padlock = rarm1 + c.spacerleft + c.anchor + geneInfo[gRegions[refIndex][i].lid].barcode + c.spacerright + rarm2;
	    double padlockGC = gccontent(padlock);
	    if ((padlockGC < minGC) || (padlockGC > maxGC)) continue;
	    
	    // Barcode
	    std::string bartmp = geneInfo[gRegions[refIndex][i].lid].barcode;
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
	    int32_t startpos = gRegions[refIndex][i].start + k + 1;
	    if (gRegions[refIndex][i].strand == '-') startpos = gRegions[refIndex][i].end - k - targetlen + 2;
	    ofile << geneInfo[gRegions[refIndex][i].lid].id << '\t';
	    ofile << geneInfo[gRegions[refIndex][i].lid].symbol << '\t';
	    ofile << geneInfo[gRegions[refIndex][i].lid].code << '\t';	    
	    ofile << c.chrname[refIndex] << ':' << startpos << '\t';
	    ofile << "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=" << c.chrname[refIndex] << ":" << startpos << "-" << startpos + targetlen - 1 << '\t'; 
	    ofile << gRegions[refIndex][i].strand << '\t';
	    ofile << c.chrname[refIndex] << ':' << gRegions[refIndex][i].start + 1 << '-' << gRegions[refIndex][i].end + 1 << '\t';
	    ofile << arm1 << '-' << arm2 << '\t';
	    ofile << c.spacerleft << '\t' << c.anchor << '\t';
	    ofile << geneInfo[gRegions[refIndex][i].lid].barcode << '\t';
	    ofile << c.spacerright << '\t';
	    ofile << padlock << '\t';
	    ofile << arm1TM << '\t' << arm2TM << '\t' << barTM << '\t' << probeTM << '\t';
	    ofile << arm1GC << '\t' << arm2GC << '\t' << barGC << '\t' << probeGC << std::endl;

	    if (c.json) {
	      if (!firstRec) rcfile << ',';
	      else firstRec = false;
	      rcfile << "[";
	      rcfile << "\"" << geneInfo[gRegions[refIndex][i].lid].id << "\", ";
	      rcfile << "\"" << geneInfo[gRegions[refIndex][i].lid].symbol << "\", ";
	      rcfile << "\"" << geneInfo[gRegions[refIndex][i].lid].code << "\", ";
	      rcfile << "\"" << c.chrname[refIndex] << ':' << startpos << "\", ";
	      rcfile << "\"" << "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=" << c.chrname[refIndex] << ":" << startpos << "-" << startpos + targetlen - 1 << "\", ";
	      rcfile << "\"" << gRegions[refIndex][i].strand << "\", ";
	      rcfile << "\"" << c.chrname[refIndex] << ':' << gRegions[refIndex][i].start + 1 << '-' << gRegions[refIndex][i].end + 1 << "\", ";
	      rcfile << "\"" << arm1 << '-' << arm2 << "\", ";
	      rcfile << "\"" << c.spacerleft << "\", ";
	      rcfile << "\"" << c.anchor << "\", ";
	      rcfile << "\"" << geneInfo[gRegions[refIndex][i].lid].barcode << "\", ";
	      rcfile << "\"" << c.spacerright << "\", ";
	      rcfile << "\"" << padlock << "\", ";
	      rcfile << "\"" << arm1TM << "\", ";
	      rcfile << "\"" << arm2TM << "\", ";
	      rcfile << "\"" << barTM << "\", ";
	      rcfile << "\"" << probeTM << "\", ";
	      rcfile << "\"" << arm1GC << "\", ";
	      rcfile << "\"" << arm2GC << "\", ";
	      rcfile << "\"" << barGC << "\", ";
	      rcfile << "\"" << probeGC << "\"";
	      rcfile << ']';
	    }
	    
	    // Increase k for non-overlapping probes
	    if (!c.overlapping) k += targetlen - 1;
	  }
	}
	free(seq);
      }
    }
    fai_destroy(fai);
    ofile.close();
    if (c.json) {
      rcfile << "]}}";
      rcfile.pop();
      rcfile.pop();
    }
    
    // Clean-up
    primer3thal::destroy_thal_structures();
    

#ifdef PROFILE
    ProfilerStop();
#endif

    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
    return 0;
  }

  inline int
  errorMessage(PadlockConfig const& c, std::string const& errmsg) {
    if (c.json) {
      boost::iostreams::filtering_ostream rcfile;
      rcfile.push(boost::iostreams::gzip_compressor());
      rcfile.push(boost::iostreams::file_sink(c.outfile.c_str(), std::ios_base::out | std::ios_base::binary));
      rcfile << "{";
      rcfile << "\"errors\": [";
      nlohmann::json err;
      err["type"] = "error";
      err["title"] = errmsg;
      rcfile << err.dump();
      rcfile << "]}";
      rcfile.pop();
      rcfile.pop();
    } else {
      std::cerr << errmsg << std::endl;
    }
    return 1;
  }
  

  int padlock(int argc, char** argv) {
    PadlockConfig c;
    
    // CMD Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
      ("gtf,t", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3 file")
      ("config,i", boost::program_options::value<boost::filesystem::path>(&c.primer3Config)->default_value("./src/primer3_config/"), "primer3 config directory")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.tsv"), "output file")
      ("json,j", boost::program_options::value<boost::filesystem::path>(&c.jsonfile), "gzipped JSON file [optional]")
      ("hamming,n", "use hamming neighborhood instead of edit distance")
      ;


    boost::program_options::options_description padopt("Padlock options");
    padopt.add_options()
      ("anchor,a", boost::program_options::value<std::string>(&c.anchor)->default_value("TGCGTCTATTTAGTGGAGCC"), "anchor sequence")
      ("spacerleft,l", boost::program_options::value<std::string>(&c.spacerleft)->default_value("TCCTC"), "spacer left")
      ("spacerright,r", boost::program_options::value<std::string>(&c.spacerright)->default_value("TCTTT"), "spacer right")
      ("barcodes,b", boost::program_options::value<boost::filesystem::path>(&c.barcodes), "FASTA barcode file")
      ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
      ("armlen,m", boost::program_options::value<uint32_t>(&c.armlen)->default_value(20), "probe arm length")
      ("probe,p", "apply distance to entire probe, i.e., only one arm needs to be unique")
      ("overlapping,v", "allow overlapping probes")
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
    cmdline_options.add(generic).add(padopt).add(tmcalc).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(padopt).add(tmcalc);
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
    if (vm.count("json")) c.json = true;
    else c.json = false;
    if (vm.count("probe")) c.armMode = false;
    else c.armMode = true;
    if (vm.count("overlapping")) c.overlapping = true;
    else c.overlapping = false;

    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      return errorMessage(c, "Error: Genome does not exist!");
    }

    // Check GTF file
    if (!(boost::filesystem::exists(c.gtfFile) && boost::filesystem::is_regular_file(c.gtfFile) && boost::filesystem::file_size(c.gtfFile))) {
      return errorMessage(c, "Error: GTF file does not exist!");
    }

    // Check barcodes file
    if (!(boost::filesystem::exists(c.barcodes) && boost::filesystem::is_regular_file(c.barcodes) && boost::filesystem::file_size(c.barcodes))) {
      return errorMessage(c, "Error: Barcode FASTA file does not exist!");
    }

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

    // Check gene list
    c.computeAll = false;
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      if (c.infile.string() == "all") c.computeAll = true;
      else c.geneset.insert(c.infile.string());
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
	    c.geneset.insert(geneName);
	  }
	}
	geneFile.close();
      }
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "dicey ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    return runPadlock(c);
  }
    
}

#endif
