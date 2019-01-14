/*
============================================================================
Silica: In-silico PCR
============================================================================
Copyright (C) 2019 Tobias Rausch

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

#ifndef ALIGN_H
#define ALIGN_H

#include <iostream>

namespace dicey
{

  template<typename TScoreValue>
  struct DnaScore {
    typedef TScoreValue TValue;

    TScoreValue match;
    TScoreValue mismatch;
    TScoreValue go;
    TScoreValue ge;
    TScoreValue inf;

    DnaScore() {
      match = 5;
      mismatch = -4;
      go = -10;
      ge = -1;
      inf = 1000000;
    }

    DnaScore(TScoreValue m, TScoreValue mm, TScoreValue gapopen, TScoreValue gapextension) : match(m), mismatch(mm), go(gapopen), ge(gapextension) {
      inf = 1000000;
    }
  };

  

  // Configure the DP matrix
  template<bool THorizontal = false, bool TVertical = false>
    class AlignConfig;

  template<>
    class AlignConfig<false, false> {};

  template<>
    class AlignConfig<false, true> {};
  
  template<>
    class AlignConfig<true, false> {};
  
  template<>
    class AlignConfig<true, true> {};

  template<bool THorizontal, typename TPos1, typename TPos2, typename TCost>
    inline TCost
    _verticalGap(AlignConfig<THorizontal, false> const&, TPos1 const, TPos2 const, TCost const cost) 
  {
    return cost;
  }
  
  template<bool THorizontal, typename TPos1, typename TPos2, typename TCost>
    inline TCost
    _verticalGap(AlignConfig<THorizontal, true> const&, TPos1 const i, TPos2 const iend, TCost const cost) 
  {
    if ((i == (TPos1) 0) || (i == (TPos1) iend)) return 0;
    else return cost;
  }

  template<bool TVertical, typename TPos1, typename TPos2, typename TCost>
    inline TCost
    _horizontalGap(AlignConfig<false, TVertical> const&, TPos1 const, TPos2 const, TCost const cost)
  {
    return cost;
  }
  
  template<bool TVertical, typename TPos1, typename TPos2, typename TCost>
    inline TCost
    _horizontalGap(AlignConfig<true, TVertical> const&, TPos1 const i, TPos2 const iend, TCost const cost)
  {
    if ((i == (TPos1) 0) || (i == (TPos1) iend)) return 0;
    else return cost;
  }

  template<typename TChar, typename TDimension>
  inline std::size_t
  _size(boost::multi_array<TChar, 2> const& a, TDimension const i) {
    return a.shape()[i];
  }

  template<typename TDimension>
  inline std::size_t
  _size(std::string const& s, TDimension const i) {
    if (i) return s.size();
    return 1;
  }


  template<typename TAIndex, typename TScore>
  inline int
  _score(std::string const& s1, std::string const& s2, TAIndex row, TAIndex col, TScore const& sc)
  {
    return (s1[row] == s2[col] ? sc.match : sc.mismatch );
  }

  template<typename TTrace, typename TAlign>
  inline void
  _createAlignment(TTrace const& trace, std::string const& s1, std::string const& s2, TAlign& align)
  {
    align.resize(boost::extents[2][trace.size()]);
    std::size_t row = 0;
    std::size_t col = 0;
    std::size_t ai = 0;
    for(typename TTrace::const_reverse_iterator itT = trace.rbegin(); itT != trace.rend(); ++itT, ++ai) {
      if (*itT == 's') {
	align[0][ai] = s1[row++];
	align[1][ai] = s2[col++];
      } else if (*itT =='h') {
	align[0][ai] = '-';
	align[1][ai] = s2[col++];
      } else {
	align[0][ai] = s1[row++];
	align[1][ai] = '-';
      }
    }
  }
}

#endif
