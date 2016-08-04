#ifndef _FASTX_IO_NUCFREQUENCY_H_
#define _FASTX_IO_NUCFREQUENCY_H_

#include <string>
#include <fstream>
#include <map>

namespace FASTX {

  // Forward declarations
  class Record;
  
  /** 
   * @brief Helper class to count nucleotide frequencies
   *
   * The NucFrequency class here is a helper class to compute
   * nucleotide frequency tables. Currently, it is implemented
   * with `std::map` as a red-black BST. 
   */
  class NucFrequency {
  public:
    /**
     * @brief Printing method (prints: Nucleotide<TAB>Count)
     *
     * @param outstream A reference to an output stream
     * @param nuc A NucFrequency object
     */
    friend std::ostream& operator<<(std::ostream& outstream, const NucFrequency nuc);
    
    /**
     * @brief Add all nucleotides in a string to the table
     *
     * @param str A reference to a string
     */
    void add(const std::string& str);

    /**
     * @brief Add all nucleotides of a `Record` to the table
     *
     * @param str A reference to a `Record` object
     */
    void add(const Record& str);
    
    /**
     * @brief Extract data for a letter
     *
     * @param nuc The letter of interest
     * @return Counts of the letter of interest
     */
    length_t operator[](char nuc) const;

    /**
     * @brief Get a vector of the letters in the count table
     *
     * @return A vector of the letters present in the count table
     */
    std::vector<char> letters(void) const;
    
    /**
     * @brief NULL constructor
     */
    NucFrequency(void) {};
    
  private:
    std::map<char, length_t> _freq_table;
  };

  /**
   * @brief Printing method (prints: Nucleotide<TAB>Count)
   *
   * @param outstream A reference to an output stream
   * @param nuc A NucFrequency object
   */
  std::ostream& operator<<(std::ostream& outstream, const NucFrequency nuc);

  /**
   * @brief Examples
   *
   * @example count_freq.cpp
   */
}
#endif
