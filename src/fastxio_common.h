/**
 * Common data variables and macros for fastx_io
 *
 * This file contains global dictionaries for nuleotide to
 * amino acid conversation, reverse complementations etc.
 */

#ifndef _FASTX_IO_COMMON_H_
#define _FASTX_IO_COMMON_H_

#include <string>
#include <set>
#include <map>
#include <vector>

#define FASTX_SW_MATCH 1
#define FASTX_SW_MISMATCH 2
#define FASTX_SW_GAP_EXTENSION 2
#define FASTX_SW_GAP_OPEN 5

#define NULL_SEQTYPE 0

#define FASTQ_TYPE 1
#define FASTA_TYPE 2

#define DNA_SEQTYPE 4
#define RNA_SEQTYPE 8
#define AA_SEQTYPE 16

#ifdef __GNU__
#ifndef PARALLEL_SORT
#include <parallel/algorithm>
#define SORT __gnu_parallel::sort
#define PARALLEL_SORT
#endif
#else
#ifndef SINGLE_SORT
#include <algorithm>
#define SORT std::sort
#define SINGLE_SORT
#endif
#endif

/**
* Root namespace
*/
namespace FASTX {

/**
 * @brief Type to store length information of sequences
 */
typedef unsigned long length_t;

/**
 * @brief Type to store scores in Smith Waterman matrix
 */
typedef long int score_t;

template<typename T>
double lexical_double(T arg)
{
  return std::stod(std::to_string(arg));
}

/**
* @brief This global struct holds all translation/complementation tables
*/
struct GData
{
  /**
   * @brief All nucleotides (IUPAC notation)
   */
  static const std::set<char> nuc_alphabet;

  /**
   * @brief All amino acid codes (IUPAC notation)
   */
  static const std::set<char> aa_alphabet;

  /**
   * @brief Reverse complementation table
   */
  static const std::map<char, char> rc;

  /**
   * @brief Disambiguation for ambiguous IUPAC DNA/RNA codes
   */
  static const std::map<char, std::vector<char> > enum_iupac_dna;

  /**
   * @brief Disambiguation for ambiguous IUPAC RNA codes
   */
  static const std::map<char, std::vector<char> > enum_iupac_rna;

  /**
   * @brief Triplett nucleotide to AA translation (DNA)
   */
  static const std::map<std::string, char> codon_to_protein_dna;

  /**
   * @brief Triplett nucleotide to AA translation (RNA)
   */
  static const std::map<std::string, char> codon_to_protein_rna;
};//struct GData

} // namespace FASTX
#endif
