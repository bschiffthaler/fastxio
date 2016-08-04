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

#define FASTQ_TYPE 1
#define FASTA_TYPE 2

#define DNA_SEQTYPE 4
#define RNA_SEQTYPE 8
#define AA_SEQTYPE 16

namespace FASTX {
  
  /** 
   * @brief Type to store length information of sequences 
   */
  typedef unsigned long length_t;
  
  /** 
   * @brief Type to store scores in Smith Waterman matrix 
   */
  typedef long int score_t; 


    /**
     * @brief All nucleotides (IUPAC notation)
     */
    static const std::set<char> G_nuc_alphabet =  {
      // Nucleotide alphabet
      'A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n',
      'U', 'u', 'R', 'r', 'Y', 'y', 'K', 'k', 'M', 'm',
      'S', 's', 'W', 'w', 'B', 'b', 'D', 'd', 'H', 'h',
      'V', 'v', '-'};

    /**
     * @brief All amino acid codes (IUPAC notation)
     */
    static const std::set<char> G_aa_alphabet =  {
      // Protein alphabet
      'A', 'a', 'B', 'b', 'C', 'c', 'D', 'd', 'E', 'e',
      'F', 'f', 'G', 'g', 'H', 'h', 'I', 'i', 'J', 'j',
      'K', 'k', 'L', 'l', 'M', 'm', 'N', 'n', 'O', 'o',
      'P', 'p', 'Q', 'q', 'R', 'r', 'S', 's', 'T', 't',
      'U', 'u', 'V', 'v', 'W', 'w', 'Y', 'y', 'Z', 'z',
      'X', 'x', '*', '-', '.'};

    /**
     * @brief Reverse complementation table
     */
    static const std::map<char, char> G_rc =  {
      {'A', 'T'}, {'a', 't'}, {'C', 'G'}, {'c', 'g'},
      {'G', 'C'}, {'g', 'c'}, {'T', 'A'}, {'t', 'a'},
      {'N', 'N'}, {'n', 'n'}, {'U', 'A'}, {'u', 'a'},
      //Special IUPAC cases
      {'R', 'Y'}, {'r', 'y'}, {'Y', 'R'}, {'y', 'r'},
      {'K', 'M'}, {'k', 'm'}, {'M', 'K'}, {'m', 'k'},
      {'S', 'S'}, {'s', 's'}, {'W', 'W'}, {'w', 'w'},
      {'B', 'V'}, {'b', 'v'}, {'V', 'B'}, {'v', 'b'},
      {'D', 'H'}, {'d', 'h'}, {'H', 'D'}, {'h', 'd'},
      {'-', '-'}
    };

    /**
     * @brief Disambiguation for ambiguous IUPAC DNA/RNA codes
     */
    static const std::map<char, std::vector<char> > G_enum_iupac_dna =  {
      {'R', {'A', 'G'}}, {'Y', {'C', 'T'}}, {'K', {'G', 'T'}},
      {'M', {'A', 'C'}}, {'S', {'C', 'G'}}, {'W', {'A', 'T'}},
      {'B', {'C', 'G', 'T'}}, {'D', {'A', 'G', 'T'}},
      {'H', {'A', 'C', 'T'}}, {'V', {'A', 'C', 'G'}},
      {'N', {'A', 'C', 'T', 'G'}}
    };

    /**
     * @brief Disambiguation for ambiguous IUPAC RNA codes
     */
    static const std::map<char, std::vector<char> > G_enum_iupac_rna =  {
      {'R', {'A', 'G'}}, {'Y', {'C', 'U'}}, {'K', {'G', 'U'}},
      {'M', {'A', 'C'}}, {'S', {'C', 'G'}}, {'W', {'A', 'U'}},
      {'B', {'C', 'G', 'U'}}, {'D', {'A', 'G', 'U'}},
      {'H', {'A', 'C', 'U'}}, {'V', {'A', 'C', 'G'}},
      {'N', {'A', 'C', 'U', 'G'}}
    };

    /**
     * @brief Triplett nucleotide to AA translation 
     */
    static const std::map<std::string, char> G_codon_to_protein =  {
      {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
      {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
      {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
      {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
      {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
      {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
      {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
      {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
      {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '.'}, {"TAG", '.'},
      {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
      {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'N'},
      {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
      {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '.'}, {"TGG", 'W'},
      {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
      {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
      {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'},
      {"UUU", 'F'}, {"UUC", 'F'}, {"UUA", 'L'}, {"UUG", 'L'},
      {"CUU", 'L'}, {"CUC", 'L'}, {"CUA", 'L'}, {"CUG", 'L'},
      {"AUU", 'I'}, {"AUC", 'I'}, {"AUA", 'I'}, {"AUG", 'M'},
      {"GUU", 'V'}, {"GUC", 'V'}, {"GUA", 'V'}, {"GUG", 'V'},
      {"UCU", 'S'}, {"UCC", 'S'}, {"UCA", 'S'}, {"UCG", 'S'},
      {"CCU", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
      {"ACU", 'U'}, {"ACC", 'U'}, {"ACA", 'U'}, {"ACG", 'U'},
      {"GCU", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
      {"UAU", 'Y'}, {"UAC", 'Y'}, {"UAA", '.'}, {"UAG", '.'},
      {"CAU", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
      {"AAU", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'N'},
      {"GAU", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
      {"UGU", 'C'}, {"UGC", 'C'}, {"UGA", '.'}, {"UGG", 'W'},
      {"CGU", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
      {"AGU", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
      {"GGU", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };
}
#endif
