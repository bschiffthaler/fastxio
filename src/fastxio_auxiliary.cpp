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
#include <fastxio_auxiliary.h>
#include <fastxio_common.h>
#include <fastxio_record.h>
#include <fastxio_reader.h>

namespace FASTX {

extern GData global; /**< Global variable of translation tables */

// Get iterator to the next whitespace character
std::basic_string<char>::const_iterator get_whitespace(std::string& s)
{
  for (std::basic_string<char>::const_iterator it = s.begin(); it != s.end(); it++)
  {
    int ascii = (int) * it;
    if (ascii == 9 || ascii == 32 || ascii == 10)
    {
      return it;
    }
  }
  return s.end();
}

// Autodetect phred encoding -->DEPRECATE in favor of char * method!
unsigned short scan_phred(std::istream& instream)
{
  std::string line;
  unsigned short NR = 0;
  while (std::getline(instream, line))
  {
    ++NR;
    if (NR % 4 == 0)
    {
      for (auto it = line.begin(); it != line.end(); it++)
      {
        if (*it < 59)
          return 33;
        if (*it > 73)
          return 64;
      }
    }
  }
  return 0;
}

unsigned short scan_phred(const char * infile)
{
  Reader r(infile, DNA_SEQTYPE);
  while (r.peek() != EOF)
  {
    Record x = r.next();
#ifndef NO_ERROR_CHEKING
    if (x.get_type() & FASTA_TYPE)
      throw std::runtime_error("FASTA files do not have quality values");
#endif
    std::string qual = x.get_qual();
    for (char c : qual)
    {
      if (c < 59) return 33;
      if (c > 73) return 64;
    }
  }
  return 0;
}

// Check magic numbers to see if file is compressed

bool is_gzip(const char * input)
{
  std::ifstream ipt(input, std::ios::in | std::ios::binary);
  char byte1 = ipt.get();
  char byte2 = ipt.get();
  char byte3 = ipt.get();
  bool gzip = ((byte1 == '\x1F') &&
               (byte2 == '\x8B') &&
               (byte3 == '\x08')) ? true : false;
  ipt.close();
  return gzip;
}

bool is_bzip2(const char * input)
{
  std::ifstream ipt(input, std::ios::in | std::ios::binary);
  char byte1 = ipt.get();
  char byte2 = ipt.get();
  char byte3 = ipt.get();
  bool bzip2 = ((byte1 == '\x42') &&
                (byte2 == '\x5a') &&
                (byte3 == '\x68')) ? true : false;
  ipt.close();
  return bzip2;
}


// Test if a character is allowed sequence (ACTGN)
bool is_sequence_char(char test, char seqtype = DNA_SEQTYPE)
{
  if ((seqtype & DNA_SEQTYPE) || (seqtype & RNA_SEQTYPE) )
  {
    if (global.nuc_alphabet.find(test) != global.nuc_alphabet.end())
      return true;
  }
  if (seqtype & AA_SEQTYPE)
  {
    if (global.aa_alphabet.find(test) != global.aa_alphabet.end())
      return true;
  }
  return false;
}

void recursive_iupac_enum(std::set<Record>& set, Record& rec,
                          const std::map<char, std::vector<char> >& translation_table,
                          std::set<char>& unambiguous_nuc)
{
  std::cerr << "call\n";
  bool all_unambiguous = true;
  size_t position = 0;
  for (char c : rec.get_seq())
  {
#ifndef NO_ERROR_CHECKING
    c = toupper(c);
    if (rec.get_type() & AA_SEQTYPE)
      throw std::runtime_error("Cannot enumerate amino acid record");
#endif
    if (unambiguous_nuc.find(c) == unambiguous_nuc.end())
    {
      all_unambiguous = false;
      std::string seq_copy = rec.get_seq();
      std::string id = rec.get_id();
      std::string qual;
      if (rec.get_type() & FASTQ_TYPE)
        qual = rec.get_qual();
      for (char r : translation_table.at(c))
      {
        std::string tmp_seq(seq_copy);
        tmp_seq[position] = r;
        std::string tmp_id = id + "_" + std::to_string(position + 1) + r;
        if (rec.get_type() & FASTQ_TYPE)
        {
          Record R(tmp_seq, tmp_id, qual,
                   rec.get_type() & DNA_SEQTYPE ? DNA_SEQTYPE : RNA_SEQTYPE);
          recursive_iupac_enum(set, R, translation_table, unambiguous_nuc);
        }
        else
        {
          Record R(tmp_seq, tmp_id,
                   rec.get_type() & DNA_SEQTYPE ? DNA_SEQTYPE : RNA_SEQTYPE);
          recursive_iupac_enum(set, R, translation_table, unambiguous_nuc);
        }
      }
      position++; break;
    }
    else
    {
      position++;
    }
  }
  if (all_unambiguous)
    set.insert(rec);
}
};
