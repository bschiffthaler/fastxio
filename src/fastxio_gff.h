#ifndef _FASTX_IO_GFF_H_
#define _FASTX_IO_GFF_H_

#include <string>
#include <fstream>
#include <memory>
#include <fastxio_common.h>
#include <unordered_map>
#include <vector>

namespace FASTX {

class GffRecord
{
public:
  GffRecord(std::string & line);
  std::string seqid = "";
  std::string source = "";
  std::string type = "";
  uint64_t start = 0;
  uint64_t end = 0;
  double score = 0.0;
  char strand = '.';
  char phase = '.';
  std::unordered_map<std::string, std::string> attributes;
  bool operator<(GffRecord const & other) const {
    if (seqid == other.seqid) {
      if (start == other.start) {
        return end < other.end;
      }
      return start < other.start;
    }
    return seqid < other.seqid;
  }
  bool operator==(GffRecord const & other) const {
    return seqid == other.seqid && start == other.start && end == other.end;
  }
  bool operator>(GffRecord const & other) {
    return ! ((*this) == other || (*this) < other);
  }
  friend std::ostream& operator<<(std::ostream& outstream, const GffRecord rhs);
};

class SequenceRegion
{
public:
  // Init
  SequenceRegion(std::string & xchr, uint64_t && xstart, uint64_t && xend) :
    chr(xchr), end(xend), start(xstart) {}

  std::string chr;
  uint64_t start;
  uint64_t end;
};

class Gff
{
public:
  Gff(std::string const & in_path);
  void read();
  std::vector<GffRecord> & records() {return _records;}
  SequenceRegion const & get_seqregion(std::string const & name) const {
    return _seqregions.at(name);
  }
private:
  std::unique_ptr<std::istream> in_handle;
  std::unordered_map<std::string, SequenceRegion> _seqregions;
  std::vector<GffRecord> _records;
  std::string _version;
};

} //Namespace FASTX

#endif