#ifndef _FASTX_IO_RECORD_H_
#define _FASTX_IO_RECORD_H_

#include <string>
#include <memory>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <matrix.h>
#include <fastxio_common.h>
#include <fastxio_nuc_frequency.h>

namespace FASTX {
  
  /** 
   * @brief General purpose FASTA and FASTQ input class.
   *
   * This library supports simple I/O, parsing and common operations
   * on FASTQ and FASTA files (abbreviated as FASTX). It's main 
   * emphases lie on *speed* and *correctness*. All IUPAC characters
   * are supported, as are upper- and lowercase notation.
   */
  class Record {
  public:
    /**
     * @brief Get the sequence.
     *
     * @return The sequence of the object.
     */
    std::string get_seq(void) const { return _seq; }

    /**
     * @brief Get the quality as ACII ancoded characters.
     *
     * @return The quality values.
     */
    std::string get_qual(void) const { return _qual; }

    /**
     * @brief Get the ID.
     * 
     * @return The ID.
     */
    std::string get_id(void) const { return _id; }

    /**
     * @brief Get the type of the record.
     *
     * The return of this function is a byte (char), which
     * contains the bit-encoded type of the record. The 
     * bits are as follows:
     * Bit | Meaning |
     * ---:|:--------|
     * 00001 | The record is FASTQ formatted |
     * 00010 | The record is FASTA formatted |
     * 00100 | The record is a DNA sequence |
     * 01000 | The record is an RNA sequence |
     * 10000 | The record is an amino acid sequence |
     * 
     * As a default, records are initialized as 00101, i.e.
     * FASTQ/DNA.
     *
     * @return The bit encoded type of the record
     */
    char get_type(void) const { return _type; }

    /**
     * @brief Get a shared pointer to sequence.
     *
     * @return A shared pointer to the sequence of the object.
     */    
    std::shared_ptr<const std::string> get_seq_ptr(void) const {
      return std::make_shared<const std::string>(_seq); }

    /**
     * @brief Get a shared pointer to quality.
     *
     * @return A shared pointer to the quality string of the object.
     */
    std::shared_ptr<const std::string> get_qual_ptr(void) const {
      return std::make_shared<const std::string>(_qual); }
    
    /**
     * @brief Get the length of the record.
     *
     * @return The length of the record.
     */
    length_t size(void) const {return _seq.length(); }

    /**
     * @brief Print the record to a sink.
     *
     * The record can be printed to a sink derived from an
     * std::ostream class. FASTA/FASTQ formatting is handled
     * automatically via the type of the read.
     *
     * @param outstream A sink to print to.
     * @param rec A record object
     */
    friend std::ostream& operator<<(std::ostream& outstream, const Record& rec);

    friend class Wrap;
    friend class NucFrequency;
    
    /**
     * @brief Constructor from an istream.
     *
     * @param input An std::istream derived source to read from
     * @param seqtype The sequence type (macro), DNA_SEQTYPE, RNA_SEQTYPE, AA_SEQTYPE
     */
    Record(std::istream& input, char seqtype);

    /**
     * @brief Constructor (FASTA) from sequence and ID.
     *
     * @param fasta The sequence for the record
     * @param id The ID for the record
     * @param seqtype The sequence type (macro), DNA_SEQTYPE, RNA_SEQTYPE, AA_SEQTYPE
     */
    Record(const std::string& fasta, const std::string& id,
	   char seqtype);

    /**
     * @brief Constructor (FASTQ) from sequence, ID and qual.
     *
     * @param fastq The sequence for the record
     * @param id The ID for the record
     * @param qual The quality string
     * @param seqtype The sequence type (macro), DNA_SEQTYPE, RNA_SEQTYPE, AA_SEQTYPE
     */
    Record(const std::string& fastq, const std::string& id, 
	   const std::string& qual, char seqtype);

    /**
     * @brief Translate a record to an AA record, one ORF. 
     *
     * @param orf The open reading frame: 0 for the first base, 1 or 2 for the next
     * @return A translated (FASTA) record
     */
    Record translate(unsigned short orf);

    /**
     * @brief Translate a record to an AA record, all ORFs.
     *
     * See also the example `translate.cpp`.
     *
     * @return A vector of all translated ORFs.
     */
    std::vector<Record> translate(void);

    /**
     * @brief Tabulate all nucleotides and get the frequencies.
     *
     * @return A class containing nucleotide frequencies.
     */
    NucFrequency count_freq(void);

    /**
     * @brief Get a sub-sequence.
     *
     * @param start The (0-offset) start position of the subsequence.
     * @param stop The stop position.
     * @return The created subsequence.
     */
    Record subseq(length_t start, length_t stop) const;

    /**
     * @brief Reverse complement a record.
     *
     * @return A reverse complement of the record.
     */
    Record rc(void) const;

    /**
     * @brief Get all k-mers of length k.
     *
     * @param k The k-mer length.
     * @return A vector containing all k-mer records.
     */
    std::vector<Record> kmer(length_t k) const;

    /**
     * @brief Create a sliding window along the record.
     *
     * @param width The width of the sliding window.
     * @param increment The increment of the window.
     * @return A vector of the sliding windows.
     */
    std::vector<Record> window(length_t width, length_t increment) const;

    /**
     * @brief Get a numeric representation of the quality values.
     *
     * The function returns a vector of numeric values of the 
     * ASCII encoded PHRED scores. The values are still containing the
     * PHRED offset. See the `scan_phred()` function in
     * fastxio_auxiliary.
     *
     * @return A vector of numeric quality values
     */
    std::vector<unsigned short> get_numeric_qual(void) const;

    /**
     * @brief Enumerate all possible sequences from ambiguous IUPAC sequences
     *
     * See also the `enumerate_iupac` example. 
     *
     * @return A set of all possible unambiguous sequences
     */
    std::set<Record> enumerate_iupac(void);

    /**
     * @brief Validate the record.
     *
     * The sequence and quality lengths need to be the same length,
     * the quality values cannot be outside the allowed range,
     * the sequence is tested for disallowed characters.
     *
     * @return True, if the record is valid, false otherwise
     */
    bool validate(void) const;

    /**
     * @brief Comparison less than operator for `Record` classes.
     *
     * The ordering of `Record` objects is alphabetical by
     * their respective sequences.
     *
     * @param a First record.
     * @param b Second record.
     *
     * @return True if the sequence of a is alphabetically before b,
     *         false otherwise.
     */
    friend bool operator<(const Record& a, const Record& b);

    /**
     * @brief Comparison greater than operator for `Record` classes.
     *
     * The ordering of `Record` objects is alphabetical by
     * their respective sequences.
     *
     * @param a First record.
     * @param b Second record.
     *
     * @return True if the sequence of a is alphabetically after b,
     *         false otherwise.
     */
    friend bool operator>(const Record& a, const Record& b);

    /**
     * @brief Comparison equal operator for `Record` classes.
     *
     * The ordering of `Record` objects is alphabetical by
     * their respective sequences.
     *
     * @param a First record.
     * @param b Second record.
     *
     * @return True if the sequence of a is identical to b,
     *         false otherwise.
     */
    friend bool operator==(const Record& a, const Record& b);

    /**
     * @brief Comparison not equal operator for `Record` classes.
     *
     * The ordering of `Record` objects is alphabetical by
     * their respective sequences.
     *
     * @param a First record.
     * @param b Second record.
     *
     * @return True if the sequence of a is different than b,
     *         false otherwise.
     */
    friend bool operator!=(const Record& a, const Record& b);

    /**
     * @brief Comparison greater than or equal operator for `Record` classes.
     *
     * The ordering of `Record` objects is alphabetical by
     * their respective sequences.
     *
     * @param a First record.
     * @param b Second record.
     *
     * @return True if the sequence of a is alphabetically after b,
     *         or identical to b, false otherwise
     */
    friend bool operator>=(const Record& a, const Record& b);

    /**
     * @brief Comparison less than or equal operator for `Record` classes.
     *
     * The ordering of `Record` objects is alphabetical by
     * their respective sequences.
     *
     * @param a First record.
     * @param b Second record.
     *
     * @return True if the sequence of a is alphabetically before b,
     *         or identical to b, false otherwise
     */
    friend bool operator<=(const Record& a, const Record& b);
  private:
    std::string _seq;
    std::string _qual;
    std::string _id;
    char _type;
    /**
     * @brief Examples
     * @example enumerate_iupac.cpp
     * @example reverse_complement.cpp
     * @example generate_kmers.cpp
     * @example generate_window.cpp
     * @example translate.cpp
     * @example numeric_qual.cpp
     */
  };

  /**
   * @brief Comparison less than operator for `Record` classes.
   */
  inline bool operator<(const Record& a, const Record& b)
  {
    if(a._seq < b._seq)
      return true;
    else if(a._seq > b._seq)
      return false;
    else
      return a._id < b._id;
  }

  /**
   * @brief Comparison greater than operator for `Record` classes.
   */
  inline bool operator>(const Record& a, const Record& b)
  {
    if(a._seq > b._seq)
      return true;
    else if(a._seq < b._seq)
      return false;
    else
      return a._id > b._id;
  }

  /**
   * @brief Comparison equal operator for `Record` classes.
   */
  inline bool operator==(const Record& a, const Record& b)
  {
    return (a._seq == b._seq && a._id == b._id);
  }

  /**
   * @brief Comparison not equal operator for `Record` classes.
   */
  inline bool operator!=(const Record& a, const Record& b)
  {
    return ! (a == b);
  }

  /**
   * @brief Comparison greater than or equal operator for `Record` classes.
   */
  inline bool operator>=(const Record& a, const Record& b)
  {
    return a._seq >= b._seq;
  }

  /**
   * @brief Comparison less than or equal operator for `Record` classes.
   */
  inline bool operator <=(const Record& a, const Record& b)
  {
    return a._seq <= b._seq;
  }

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
  class Wrap
  {
  public:
    
    /**
     * @brief Constructor from const Record
     *
     * @param rec A record object
     */
    Wrap(const Record& rec);

    /**
     * @brief Constructor from const Record with specified width
     *
     * @param rec A record object
     * @param width The column width (Default: 80)
     */
    Wrap(const Record& rec, unsigned int width);
    
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
    friend std::ostream& operator<<(std::ostream& outstream, Wrap rec);

  private:
    std::shared_ptr<const std::string> _id;
    std::shared_ptr<const std::string> _seq;
    unsigned int _width = 80;
    
    /**
     * Examples
     *  @example print_wrapped.cpp
     */
  };

  /**
   * @brief Print the record to a sink.
   *
   * The record can be printed to a sink derived from an
   * std::ostream class. FASTA/FASTQ formatting is handled
   * automatically via the type of the read.
   *
   * @param outstream A sink to print to.
   * @param rec A record object
   */
  std::ostream& operator<<(std::ostream& outstream, const Record& rec);

  /**
   * @brief Print a FASTA record wrapped to a specified column number (internal function)
   *
   * @param outstream The output stream
   * @param rec The wrapped Record object
   */
  std::ostream& operator<<(std::ostream& outstream, Wrap rec);
  
}
#endif
