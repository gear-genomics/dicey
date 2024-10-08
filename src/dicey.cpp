#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "bed.h"
#include "version.h"
#include "index.h"
#include "blacklist.h"
#include "chop.h"
#include "hunter.h"
#include "silica.h"
#include "mapbam.h"
#include "padlock.h"

using namespace dicey;

inline void
displayUsage() {
  std::cout << "Usage: dicey <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Commands:" << std::endl;
  std::cout << std::endl;
  std::cout << "    index        index FASTA reference file" << std::endl;
  std::cout << "    hunt         search DNA sequences" << std::endl;
  std::cout << "    search       in-silico PCR" << std::endl;
  std::cout << "    padlock      padlock probe design" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Mappability:" << std::endl;
  std::cout << std::endl;
  std::cout << "    chop         chop reference into paired-end / single-end reads" << std::endl;
  std::cout << "    mappability2 parse BAM from mapped chopped reads (requires chop + map before)" << std::endl;
  std::cout << "    blacklist    blacklist certain regions in mappability map" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

int main(int argc, char **argv) {
  if (argc < 2) { 
    printTitle("Dicey");
    displayUsage();
    return 0;
  }
  
  if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    std::cout << "Dicey version: v" << diceyVersionNumber << std::endl;
    std::cout << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
    std::cout << " using HTSlib: v" << hts_version() << std::endl;
    return 0;
  }
  else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    printTitle("Dicey");
    displayUsage();
    return 0;
  }
  else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
    displayWarranty();
    return 0;
  }
  else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
    bsd();
    return 0;
  }
  else if ((std::string(argv[1]) == "index")) {
    return index(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "hunt")) {
    return hunter(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "padlock")) {
    return padlock(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "search")) {
    return silica(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "mappability2")) {
    return mapbam(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "chop")) {
    return chop(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "blacklist")) {
    return blacklist(argc-1,argv+1);
  } else {
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
  }
}

