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
#include <fastxio_record.h>
#include <fastxio_reader.h>
#include <fastxio_auxiliary.h>

namespace FASTX {

// Open filepath automatically and detect compression
Reader::Reader(const char* infile, const char seqtype) :
  _seqtype(seqtype)
{
  if (is_gzip(infile))
  {
    _istream.reset(new igzstream(infile));
  }
  else if (is_bzip2(infile))
  {
    _istream.reset(new ibz2stream(infile));
  }
  else
  {
    _istream.reset(new std::ifstream(infile));
  }
#ifndef NO_ERROR_CHECKING
  if (! _istream->good())
  {
    throw std::runtime_error("Could not open file: " + std::string(infile));
  }
#endif
}

// Get next record
Record Reader::next(void)
{
  return Record(*_istream, _seqtype);
}

// Get next character
char Reader::peek(void)
{
  return _istream->peek();
}

};
