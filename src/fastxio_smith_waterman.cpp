#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#ifndef NO_ERROR_CHECKING
  #include <cerrno>
  #include <stdexcept>
#endif

#include <gzstream.h>
#include <bz2stream.h>
#include <matrix.h>
#include <set>
#include <map>
#include <fastxio_common.h>
#include <fastxio_smith_waterman.h>


namespace FASTX {

  score_t sw_max_neighbour(score_t a, score_t b, score_t c)
  {
    if(a >= b && a >= c) return a;
    if(b >= a && b >= c) return b;
    if(c >= a && c >= b) return c;
    return 0;
  }

  SW::SW(std::string s1, std::string s2)
  {
    //Shorter sequence always s2
    if(s1.size() >= s2.size())
      {
        _s1 = s1; _s2 = s2;
      }
    else
      {
        _s1 = s2; _s2 = s1;
      }

    _scores = Matrix<score_t>(_s2.size() + 1, _s1.size() + 1);

    for(unsigned int j = 0; j < _s1.size() + 1; j++)
      {
        _scores.at(0,j) = 0;
      }
    for(unsigned int i = 0; i < _s2.size() + 1; i++)
      {
        _scores.at(i,0) = 0;
      }

    SW_path max_idx;
    score_t Mmax = 0;

    bool previous_was_gap = false;

    for(unsigned int i = 1; i < _s2.size() + 1; i++)
      {
        for(unsigned int j = 1; j < _s1.size() + 1; j++)
          {
            score_t a = _s1[j - 1] == _s2[i - 1] ?
              _scores.at(i - 1, j - 1) + FASTX_SW_MATCH :
              _scores.at(i - 1, j - 1) - FASTX_SW_MISMATCH;
            score_t b;
            score_t c;
            if(previous_was_gap)
              {
                b = _scores.at(i - 1, j) - FASTX_SW_GAP_OPEN;
                c = _scores.at(i, j - 1) - FASTX_SW_GAP_OPEN;
              }
            else
              {
                b = _scores.at(i - 1, j) - FASTX_SW_GAP_EXTENSION;
                c = _scores.at(i, j - 1) - FASTX_SW_GAP_EXTENSION;
              }
            score_t s = sw_max_neighbour(a, b, c);
            if(s != a)
              {
                previous_was_gap = true;
              }
            else
              {
                previous_was_gap = false;
              }
            _scores.at(i, j) = s > 0 ? s : 0;

            if(s >= Mmax)
              {
                Mmax = s;
                max_idx.j = j;
                max_idx.i = i;
              }
          }
      }


    _path.push_back(max_idx);
    for(unsigned int i = max_idx.i + 1; i > 0; i--)
      {
        SW_path x;
        score_t a = _scores.at(max_idx.i - 1, max_idx.j - 1);
        score_t b = _scores.at(max_idx.i - 1, max_idx.j);
        score_t c = _scores.at(max_idx.i, max_idx.j - 1);
        if(a >= b && a >= c)
          {
            x.i = max_idx.i - 1; x.j = max_idx.j - 1;
          }
        else if(c >= a && c >= b)
          {
            x.i = max_idx.i; x.j = max_idx.j - 1;
          }
        else if(b >= a && b >= c)
          {
            x.i = max_idx.i - 1; x.j = max_idx.j;
          }

        _path.push_back(x);
        max_idx.i = x.i;
        max_idx.j = x.j;
      }
    _s1_stop = (_path.begin())->j;
    _s1_start = (_path.end() - 1)->j;
    _s2_stop = (_path.begin())->i;
    _s2_start = (_path.end() - 1)->i;

  }

  void SW::print_matrix(std::ostream& out)
  {
    for(unsigned int i = 0; i < _s2.size() + 1; i++)
      {
        if(i == 0)
          {
            out << "\t-\t";
            for(auto c : _s1) out << c << '\t';
          }
        out << '\n';
        for(unsigned int j = 0; j < _s1.size() + 1; j++)
          {
            if(j == 0)
              {
                if(i == 0)
                  {
                    out << "-\t";
                  }
                else
                  {
                    out << _s2[i - 1] << '\t';
                  }
              }
            out << _scores.at(i, j) << '\t';
          }
        out << '\n';
      }
  }

  void SW::print_path(std::ostream& out)
  {
    for(auto it = _path.begin(); it != _path.end(); it++)
      {
        out << '(' << it->i << ',' << it->j << ')';
        if(it != (_path.end() - 1)) out << "->";
      }
    out << '\n';
  }


  std::ostream& operator<<(std::ostream& out, const SW sw)
  {
    std::string x;
    std::string y;
    SW_path current = sw._path[0];
    for(auto it = (sw._path.begin() + 1); it != sw._path.end(); it++)
      {
        if(current.j - it->j == 1 && current.i - it->i == 1)
          {
            x += sw._s1[current.j - 1]; y += sw._s2[current.i - 1];
            current.j--; current.i--;
          }
        if((current.j - it->j) == 1 && (current.i - it->i) == 0)
          {
            x += sw._s1[current.j - 1]; y += '-';
            current.j--;
          }
        if((current.j - it->j) == 0 && (current.i - it->i) == 1)
          {
            x += '-'; y += sw._s2[current.i - 1];
            current.i--;
          }
      }
    out << sw._s1_start << '\t';
    for(auto it = x.rbegin(); it != x.rend(); it++)
      out << *it;
    out << '\t' << sw._s1_stop << '\n' << sw._s2_start << '\t';
    for(auto it = y.rbegin(); it != y.rend(); it++)
      out << *it;
    out << '\t' <<  sw._s2_stop << '\n';
    return out;
  }  

};
