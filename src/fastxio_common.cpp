#include <string>
#include <set>
#include <map>
#include <vector>
#include <fastxio_common.h>


namespace FASTX {


const std::set<char> GData::nuc_alphabet =  {
  'A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'N', 'n',
  'U', 'u', 'R', 'r', 'Y', 'y', 'K', 'k', 'M', 'm',
  'S', 's', 'W', 'w', 'B', 'b', 'D', 'd', 'H', 'h',
  'V', 'v', '-'
};

const std::set<char> GData::aa_alphabet =  {
  'A', 'a', 'B', 'b', 'C', 'c', 'D', 'd', 'E', 'e',
  'F', 'f', 'G', 'g', 'H', 'h', 'I', 'i', 'J', 'j',
  'K', 'k', 'L', 'l', 'M', 'm', 'N', 'n', 'O', 'o',
  'P', 'p', 'Q', 'q', 'R', 'r', 'S', 's', 'T', 't',
  'U', 'u', 'V', 'v', 'W', 'w', 'Y', 'y', 'Z', 'z',
  'X', 'x', '*', '-', '.'
};

const std::map<char, char> GData::rc =  {
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

const std::map<char, std::vector<char> > GData::enum_iupac_dna =  {
  {'R', {'A', 'G'}}, {'Y', {'C', 'T'}}, {'K', {'G', 'T'}},
  {'M', {'A', 'C'}}, {'S', {'C', 'G'}}, {'W', {'A', 'T'}},
  {'B', {'C', 'G', 'T'}}, {'D', {'A', 'G', 'T'}},
  {'H', {'A', 'C', 'T'}}, {'V', {'A', 'C', 'G'}},
  {'N', {'A', 'C', 'T', 'G'}}
};

const std::map<char, std::vector<char> > GData::enum_iupac_rna =  {
  {'R', {'A', 'G'}}, {'Y', {'C', 'U'}}, {'K', {'G', 'U'}},
  {'M', {'A', 'C'}}, {'S', {'C', 'G'}}, {'W', {'A', 'U'}},
  {'B', {'C', 'G', 'U'}}, {'D', {'A', 'G', 'U'}},
  {'H', {'A', 'C', 'U'}}, {'V', {'A', 'C', 'G'}},
  {'N', {'A', 'C', 'U', 'G'}}
};

const std::map<std::string, char> GData::codon_to_protein_dna =  {
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
  {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
};

const std::map<std::string, char> GData::codon_to_protein_rna =  {
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

GData global;

}

