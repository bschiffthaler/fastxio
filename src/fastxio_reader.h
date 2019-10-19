#ifndef _FASTX_IO_READER_H_
#define _FASTX_IO_READER_H_

#include <string>
#include <memory>
#include <fstream>
#include <fastxio_record.h>

namespace FASTX {

/**
 * @brief Class to automatically read records from files.
 *
 * This class is mainly a wrapper around the `std::istream&`
 * constructor and is able to detect compression based on the
 * magic number of the file that is passed to it.
 *
 */
class Reader {
public:
  /**
   * @brief File path constructor.
   *
   * @param file A path to a file
   * @param seqtype The sequence type: DNA_SEQTYPE, RNA_SEQTYPE, AA_SEQTYPE
   */
  Reader(const char * file, const char seqtype);

  /**
   * @brief Return next record.
   *
   * @return A record object of the next entry in the file
   */
  Record next();

  /**
   * @brief Peek the next character.
   *
   * This method is a passthrough to the `peek()` method of the
   * `std::istream` class and is mainly meant as a way to check
   * for the end of file character when reading the entire file.
   *
   * @return The next character, or EOF.
   */
  char peek(void);

  /**
   * @brief Get stream offset;
   *
   * This method passes `tellg()` to the underlying stream. Currently
   * this will only be correct for *uncompressed* streams.
   *
   * @return File offset
   */
  int tell(void) { return _istream->tellg(); }

  /**
   * @brief Seek to offset;
   *
   * This method passes `seekg()` to the underlying stream. Currently
   * this will only be correct for *uncompressed* streams.
   *
   */
  void seek(int offset) { _istream->seekg(offset); }
private:
  std::unique_ptr<std::istream> _istream;
  const char _seqtype;
};

}
#endif
