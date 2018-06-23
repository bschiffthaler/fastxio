#ifndef _FASTX_IO_NUCFREQUENCY_H_
#define _FASTX_IO_NUCFREQUENCY_H_

#include <string>
#include <fstream>
#include <map>

namespace FASTX {

  // Forward declarations
  class Record;
  class NucPercent;
  
  /** 
   * @brief Helper class to count nucleotide frequencies
   *
   * The NucFrequency class here is a helper class to compute
   * nucleotide frequency tables. Currently, it is implemented
   * with `std::map` as a red-black BST. 
   */
  class NucFrequency {
  public:
    friend class NucPercent;
    /**
     * @brief Printing method (prints: Nucleotide<TAB>Count)
     *
     * @param outstream A reference to an output stream
     * @param nuc A NucFrequency object
     */
    friend std::ostream& operator<<(std::ostream& outstream, const NucFrequency nuc);

    friend std::ostream& operator<<(std::ostream& outstream, NucPercent rhs);
    
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
     * @brief Extract data for a letter with bounds checking
     *
     * @param nuc The letter of interest
     * @return Counts of the letter of interest
     */
    length_t at(char nuc) const;

    /**
     * @brief Extract data for a letter
     *
     * @param nuc The letter of interest
     * @return Counts of the letter of interest
     */
    length_t operator[](char nuc);

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
   * @brief Helper class to wrap FASTA when printing
   *
   * Typically, FASTA records are printed to wrap lines at
   * 80 characters per column. This method does just this.
   * 
   * @warning Note that FASTQ records cannot be printed with this method
   * because they need consitently four lines per record.
   * See the example source `print_wrapped.cpp`.
   */
  class NucPercent
  {
  public:
    
    /**
     * @brief Constructor from const NucFrequency
     *
     * @param freq A NucFrequency object
     */
    NucPercent(const NucFrequency& freq): _freq(freq) {}
    
    /**
     * @brief Overload of operator() to handle formatting and passing to an std::ostream
     *
     * This method is used internally to be called from operator<<
     *
     * @param outstream The output stream
     */
    std::ostream& operator()(std::ostream& outstream) const;

    /**
     * @brief Overloaded operator<< to print to an ostream
     *
     * @param outstream The output stream
     * @param rec The wrapped Record object
     */
    friend std::ostream& operator<<(std::ostream& outstream, NucPercent rhs);
  private:
    const NucFrequency& _freq;
  };

  /**
   * @brief Printing method (prints: Nucleotide<TAB>Count)
   *
   * @param outstream A reference to an output stream
   * @param nuc A NucFrequency object
   */
  std::ostream& operator<<(std::ostream& outstream, const NucFrequency nuc);

  std::ostream& operator<<(std::ostream& outstream, NucPercent rhs);

  /**
   * @brief Examples
   *
   * @example count_freq.cpp
   */
}
#endif
