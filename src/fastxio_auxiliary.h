#ifndef _FASTX_IO_AUXILIARY_H_
#define _FASTX_IO_AUXILIARY_H_
/** 
 * @brief Auxiliary functions
 *
 * The functions here provide auxiliary features to the classes,
 * most notably detecting formats and parsing characters.
 *
*/

#include <string>
#include <fstream>
#include <set>
#include <vector>
#include <fastxio_record.h>

namespace FASTX {

  /**
   * @brief Return an iterator to the first whitespace character.
   * 
   * @param str An input string.
   * @return Iterator to the first whitespace character in the string.
   */
  std::basic_string<char>::const_iterator get_whitespace(std::string& str);

  /**
   * @brief Detect the PHRED encoding offset
   * 
   * @param instream Input stream
   * @return PHRED offset, 33 or 64
   */
  unsigned short scan_phred(std::istream& instream);

  /**
   * @brief Detect the PHRED encoding offset
   *
   * @param infile Path to the FASTQ file (can be GZip, BZip2 compressed)
   * @return PHRED offset, 33 or 64
   */
  unsigned short scan_phred(const char * infile);

  /**
   * @brief Detect if a file is GZip compressed. 
   * 
   * Specifically, this function tests
   * the first three bytes (magic number) for the zlib signature. While this
   * method is generally correct, there is a possibility for false 
   * positives/negatives.
   *
   * @param input Path to the file
   * @return True if GZip compressed, false otherwise.
   */
  bool is_gzip(const char * input);

  /**
   * @brief Detect if a file is BZip2 compressed. 
   *
   * Specifically, this function tests
   * the first three bytes (magic number) for the bzlib signature. While this
   * method is generally correct, there is a possibility for false
   * positives/negatives.
   *
   * @param input Path to the file
   * @return True if BZip2 compressed, false otherwise.
   */
  bool is_bzip2(const char * input);

  /**
   * @brief Check if a character is an allowed sequence character in DNA, RNA, or
   * amino acid sequence.
   *
   * @param test The character to be tested
   * @param seqtype The sequence type (macro) DNA_SEQTYPE, RNA_SEQTYPE or AA_SEQTYPE
   * @return True, if the character is allowed, flase otherwise.
   */
  bool is_sequence_char(const char test, char seqtype);

  /**
   * @brief Recursively enumerate all possibilities from ambiguous IUPAC DNA/RNA
   * characters. 
   *
   * This function not intended for the end user. See rather the
   * `enumerate()` method of the `Record` class.
   *
   * @param set A set containing fully enumerated sequences
   * @param rec The record to be enumerated
   * @param translation_table A table that translates ambiguous
   *        nucleotides to the possibilities
   * @param unambiguous_nuc A set that contains the unambiguous characters
   * @sa    `Record.enumerate()`
   */
  void recursive_iupac_enum(std::set<Record>& set, Record& rec,
                            const std::map<char, std::vector<char> >& translation_table,
                            std::set<char>& unambiguous_nuc);
  
}
#endif
