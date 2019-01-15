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

#ifndef NEEDLE_H
#define NEEDLE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <iostream>
#include "align.h"

namespace dicey
{

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    TMatrix mat(boost::extents[m+1][n+1]);

    // Initialization
    mat[0][0] = 0;
    for(std::size_t col = 1; col <= n; ++col) mat[0][col] = mat[0][col-1] + _horizontalGap(ac, 0, m, sc.ge);
    for(std::size_t row = 1; row <= m; ++row) mat[row][0] = mat[row-1][0] + _verticalGap(ac, 0, n, sc.ge);

    // Recursion
    for(std::size_t row = 1; row <= m; ++row)
      for(std::size_t col = 1; col <= n; ++col)
	mat[row][col] = std::max(std::max(mat[row-1][col-1] + _score(a1, a2, row-1, col-1, sc), mat[row-1][col] + _verticalGap(ac, col, n, sc.ge)), mat[row][col-1] + _horizontalGap(ac, row, m, sc.ge));

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if ((row>0) && (mat[row][col] == mat[row-1][col] + _verticalGap(ac, col, n, sc.ge))) {
	--row;
	trace.push_back('v');
      } else if ((col>0) && (mat[row][col] == mat[row][col-1] + _horizontalGap(ac, row, m, sc.ge))) {
	--col;
	trace.push_back('h');
      } else {
	--row;
	--col;
	trace.push_back('s');
      }
    }

    // Create alignment
    _createAlignment(trace, a1, a2, align);

    // Score
    return mat[m][n];
  }

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac)
  {
    DnaScore<int> dnasc;
    return needle(a1, a2, align, ac, dnasc);
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  needle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    AlignConfig<false, false> ac;
    return needle(a1, a2, align, ac);
  }


  inline void
  _debugAlignment(std::string const& pr, std::string const& genomicseq, int32_t const fwdrev) {
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    AlignConfig<true, false> semiglobal;
    DnaScore<int> sc(5, -4, -4, -4);
    int32_t score = 0;
    std::string primer(pr);
    if (fwdrev == 0) score = needle(primer, genomicseq, align, semiglobal, sc);
    else {
      reverseComplement(primer);
      score = needle(primer, genomicseq, align, semiglobal, sc);
    }
    if (fwdrev == 0) std::cerr << ">fwd";
    else std::cerr << ">rev";
    std::cerr << " AlignScore: " << score << std::endl;
    typedef typename TAlign::index TAIndex;
    for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	std::cerr << align[i][j];
      }
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }

}

#endif
