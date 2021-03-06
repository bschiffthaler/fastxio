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
#include <fastxio_nuc_frequency.h>
#include <fastxio_record.h>

namespace FASTX {

// Print nucleotide frequency
std::ostream& operator<<(std::ostream& outstream, const NucFrequency nuc)
{
  for (auto c : nuc._freq_table)
    outstream << c.first << '\t' << c.second << '\n';
  return outstream;
}

// Count nucleotides
void NucFrequency::add(const std::string& seq)
{
  for (char c : seq)
  {
    _freq_table[c]++;
  }
}

void NucFrequency::add(const Record& rec)
{
  for (char c : rec._seq)
  {
    _freq_table[c]++;
  }
}

length_t NucFrequency::at(char nuc) const
{
  return _freq_table.at(nuc);
}

length_t NucFrequency::operator[](char nuc)
{
  return _freq_table[nuc];
}

std::vector<char> NucFrequency::letters(void) const
{
  std::vector<char> ret;
  for (auto c : _freq_table)
    ret.push_back(c.first);
  return ret;
}

std::ostream& operator<<(std::ostream& outstream, NucPercent rhs)
{
  uint64_t sum = 0;
  for (auto& i : rhs._freq._freq_table)
    sum += i.second;

  double dsum = lexical_double(sum);
  for (auto& i : rhs._freq._freq_table)
  {
    outstream << i.first << '\t' << lexical_double(i.second) / dsum * 100 << '\n';
  }

  return outstream;
}

};
