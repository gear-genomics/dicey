#ifndef VERSION_H
#define VERSION_H

namespace dicey
{


  std::string diceyVersionNumber = "0.2.5";

  inline 
  void printTitle(std::string const& title) 
    {
      std::cout << "**********************************************************************" << std::endl;
      std::cout << "Program: Dicey" << std::endl;
      std::cout << "This is free software, and you are welcome to redistribute it under" << std::endl;
      std::cout << "certain conditions (GPL Version 3); for license details use '-l'." << std::endl;
      std::cout << "This program comes with ABSOLUTELY NO WARRANTY; for details use '-w'." << std::endl;
      std::cout <<  std::endl;
      std::cout <<  title << " (Version: " << diceyVersionNumber << ")" << std::endl;
      std::cout << "Contact: Tobias Rausch (rausch@embl.de)" << std::endl;
      std::cout << "**********************************************************************" << std::endl;
      std::cout << std::endl;
    }

  inline
  void displayWarranty() {
    std::cout << "Copyright (C) 2019 European Molecular Biology Laboratory (EMBL)" << std::endl;
    std::cout << std::endl;
    std::cout << "This program is free software: you can redistribute it and/or modify" << std::endl;
    std::cout << "it under the terms of the GNU General Public License as published by" << std::endl;
    std::cout << "the Free Software Foundation, either version 3 of the License, or" << std::endl;
    std::cout << "(at your option) any later version." << std::endl;
    std::cout << std::endl;
    std::cout << "This program is distributed in the hope that it will be useful," << std::endl;
    std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
    std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
    std::cout << "GNU General Public License for more details." << std::endl;
    std::cout	<< std::endl;
    std::cout << "You should have received a copy of the GNU General Public License" << std::endl;
    std::cout << "along with this program.  If not, see https://www.gnu.org/licenses/" << std::endl;
  }
    
  inline void
  bsd() {
    std::cout << "Program: Dicey" << std::endl;
    std::cout << "Version: " << diceyVersionNumber << std::endl;
    std::cout << "Copyright (C) 2019 European Molecular Biology Laboratory (EMBL) (see LICENSE file for details)." << std::endl;
    std::cout << "All rights reserved." << std::endl;
    std::cout << std::endl;
  }

}

#endif
