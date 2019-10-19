#include <string>
#include <fstream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <fastxio_auxiliary.h>
#include <fastxio_common.h>
#include <fastxio_gff.h>
#include <unordered_map>
#include <gzstream.h>
#include <bz2stream.h>
#include <str_manip.h>

namespace FASTX {

Gff::Gff(std::string const & in_path) {
  if (is_gzip(in_path.c_str()))
  {
    in_handle.reset(new igzstream(in_path.c_str()));
  }
  else if (is_bzip2(in_path.c_str()))
  {
    in_handle.reset(new ibz2stream(in_path.c_str()));
  }
  else
  {
    in_handle.reset(new std::ifstream(in_path.c_str()));
  }
  if (! in_handle->good())
  {
    throw std::runtime_error("File " + in_path + " cannot be opened.");
  }
  read();
}

void Gff::read()
{
  std::string line;
  std::string version;
  uint64_t lno = 0;
  while (std::getline(*in_handle, line))
  {
    lno++;
    if (line.front() == '#')
    {
      if (BS::str_startswith(line, "#gff-version") ||
          BS::str_startswith(line, "##gff-version"))
      {
        auto split = BS::str_split_np(line);
        _version = split[1];
      }
      else if (BS::str_startswith(line, "#sequence-region") ||
               BS::str_startswith(line, "##sequence-region"))
      {
        auto split = BS::str_split_np(line);
        _seqregions.emplace(split[1],
                            SequenceRegion(split[1],
                                           std::stoul(split[2]),
                                           std::stoul(split[3])));
      }
    }
    else
    {
      GffRecord rec(line);
      _records.push_back(rec);
    }
  }
  std::sort(_records.begin(), _records.end());
}

GffRecord::GffRecord(std::string & line)
{
  auto split = BS::str_split(line, '\t');
  if (split.size() != 9)
  {
    throw std::runtime_error("Record " + line + " does not have 9 fields");
  }
  seqid = split[0];
  source = split[1];
  type = split[2];
  start = std::stoul(split[3]);
  end = std::stoul(split[4]);
  if (split[5] == ".")
  {
    score = std::numeric_limits<double>::signaling_NaN();
  }
  else
  {
    score = std::stod(split[5]);
  }
  strand = split[6].front();
  phase = split[7].front();
  auto attrs = BS::str_split(split[8], ';');
  for (auto & attr_pair : attrs)
  {
    auto key_value = BS::str_split(attr_pair, '=');
    attributes.emplace(key_value[0], key_value[1]);
  }
}

std::ostream& operator<<(std::ostream& outstream, const GffRecord rhs)
{
  outstream <<
            rhs.seqid << '\t' <<
            rhs.source << '\t' <<
            rhs.type << '\t' <<
            rhs.start << '\t' <<
            rhs.end << '\t';
  if (std::isnan(rhs.score))
  {
    outstream << ".\t";
  }
  else
  {
    outstream << rhs.score << '\t';
  }
  outstream <<
            rhs.strand << '\t' <<
            rhs.phase << '\t';
  uint64_t kv_out = 0;
  for (auto & key_value : rhs.attributes)
  {
    outstream << key_value.first << '=' << key_value.second;
    kv_out++;
    if (kv_out < rhs.attributes.size())
    {
      outstream << ';';
    }
  }
  return outstream;
}

}
