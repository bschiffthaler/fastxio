// Copyright (2018) <BAstian Schiffthaler>
#include <omp.h>
#include <fastxio_reader.h>
#include <fastxio_record.h>
#include <fastxio_minhash.h>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

std::string to_short(const std::string& chrom) {
  std::stringstream ret;
  for (const char& c : chrom) {
    if (! std::isspace(c)) {
      ret << c;
    } else {
      break;
    }
  }
  return ret.str();
}

int main(int argc, char ** argv)
{
  try
  {
    std::string ref;  // reference sequence file
    bool one_based = false;

    po::options_description umbrella;
    po::options_description opt("Options");

    po::options_description req("Required");
    req.add_options()
    ("ref", po::value<std::string>(&ref)->required(),
     "Input reference sequence")
    ("one", po::bool_switch(&one_based),
     "Report 1-based offsets in BED output (default are 0-based offsets)")
    ;

    umbrella.add(opt).add(req);

    po::positional_options_description p;
    p.add("ref", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(umbrella).positional(p).run(), vm);

    if (vm.count("help") || argc == 1)
    {
      std::cerr << umbrella << '\n';
      return EINVAL;
    }

    po::notify(vm);

    // Create MinHash reference from first file (fasta)
    FASTX::Reader ref_reader(ref.c_str(), DNA_SEQTYPE);

    while (ref_reader.peek() != EOF) {
      FASTX::Record r = ref_reader.next();
      const std::string& s = r.get_seq();
      const std::string& chrom = r.get_id();
      uint64_t start = 0;
      for (auto it = s.cbegin(); it != s.cend();) {
        if ((*it) == 't' || (*it) == 'T') {
          it++;
          start++;
          if (it == s.cend()) {
            break;
          }
          if ((*it) == 'a' || (*it) == 'A') {
            auto left = one_based ? start : start - 1;
            auto right = one_based ? start + 1 : start;
            std::cout 
              << to_short(chrom) << '\t' 
              << left << '\t' 
              << right << '\n';
          }
        } else {
          it++;
          start++;
        }
      }
    }

  }
  catch (std::exception& e)
  {
    std::cerr << "Exception occured: " << e.what() << '\n';
    return 1;
  }

  return 0;
}
