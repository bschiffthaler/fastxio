#include <boost/program_options.hpp>
#include <iostream>
#include <fastxio_record.h>
#include <fastxio_reader.h>
#include "robin_hood.h"
#include <omp.h>

namespace po = boost::program_options;

struct param
{
  std::string genome_fasta;
  std::string kmer_fasta;
};

struct offset_t
{
  uint64_t chr;
  uint64_t pos;
};

struct k_map_t
{
  robin_hood::unordered_map<std::string, std::vector<offset_t>> map;
  std::vector<std::string> ids;
};

// Reads first record from FASTA file and gets K
uint64_t get_k(std::string const & kmer_file)
{
  FASTX::Reader R(kmer_file.c_str(), DNA_SEQTYPE);
  if (R.peek() == EOF)
  {
    throw std::runtime_error("Error opening kmer file.");
  }
  FASTX::Record r = R.next();
  return r.size();
}

k_map_t build_reference(std::string const & genome_file, uint64_t k)
{
  std::cerr << "Building index...\n";
  k_map_t ret;
  FASTX::Reader R(genome_file.c_str(), DNA_SEQTYPE);
  if (R.peek() == EOF)
  {
    throw std::runtime_error("Error opening genome file.");
  }
  uint64_t id = 0;
  uint64_t nuc = 0;
  while (R.peek() != EOF)
  {
    FASTX::Record r = R.next();
    std::cerr << "Processing: " << r.get_id() << "...\n";
    ret.ids.push_back(r.get_id());
    auto seq_ptr = r.get_seq_ptr();
    uint64_t start = 0;
    uint64_t stop = k;
    while (stop < seq_ptr->size())
    {
      std::string kmer = seq_ptr->substr(start, k);
      ret.map[kmer].push_back(offset_t({id, start + 1}));
      start++;
      stop++;
    }
    nuc += r.size();
    std::cerr << "Have " << nuc << " nucleotides\n";
    std::cerr << "Have " << ret.map.size() << " kmers\n";
    id++;
  }
  return ret;
}

void map(std::string const & kmer_file, k_map_t const & kmap)
{
  FASTX::Reader R(kmer_file.c_str(), DNA_SEQTYPE);
  std::vector<FASTX::Record> buffer;
  char const * nucs = "ACGT";

  while (R.peek() != EOF)
  {
    FASTX::Record r = R.next();
    buffer.push_back(r);
  }

  #pragma omp parallel for
  for (uint64_t i = 0; i < buffer.size(); i++)
  {
    FASTX::Record r = buffer[i];
    FASTX::Record r_rc = !r;
    auto seqptr = r.get_seq_ptr();
    auto mapptr = kmap.map.find((*seqptr));
    auto seqptr_rc = r_rc.get_seq_ptr();
    auto mapptr_rc = kmap.map.find((*seqptr_rc));

    std::vector<offset_t> hits;

    for (uint64_t n = 0; n < seqptr->size(); n++)
    {
      for (uint64_t m = 0; m < 4; m++)
      {
        std::string scopy = (*seqptr);
        scopy[n] = nucs[m];
        auto mapptr = kmap.map.find(scopy);
        if (mapptr != kmap.map.end())
        {
          for (offset_t const & hit : mapptr->second)
          {
            hits.push_back(hit);
          }
        }
      }
    }
    for (uint64_t n = 0; n < seqptr_rc->size(); n++)
    {
      for (uint64_t m = 0; m < 4; m++)
      {
        std::string scopy = (*seqptr_rc);
        scopy[n] = nucs[m];
        auto mapptr = kmap.map.find(scopy);
        if (mapptr != kmap.map.end())
        {
          for (offset_t const & hit : mapptr->second)
          {
            hits.push_back(hit);
          }
        }
      }
    }

    #pragma omp critical
    {
      if (hits.size() == 0)
      {
        std::cout <<
                  *seqptr << '\t' <<
                  "NA\t"
                  "NA\n";
      }
      else
      {
        for (offset_t const & hit : hits)
        {
          std::cout <<
                    *seqptr << '\t' <<
                    kmap.ids[hit.chr] << '\t' <<
                    hit.pos << '\n';
        }
      }
    }
  }
}

int main(int argc, char const ** argv)
{
  try
  {
    param opts;
    namespace po = boost::program_options;

    po::variables_map vm;
    po::options_description umbrella("Extract regulatory regions of genes "
                                     "upstream and downstream.");

    po::options_description opt("Optional");
    opt.add_options()
    ("help,h", "Show this help message");

    po::options_description req("Required");
    req.add_options()
    ("fasta,f", po::value<std::string>(&opts.genome_fasta)->required(),
     "Genome FASTA file")
    ("kmers,k", po::value<std::string>(&opts.kmer_fasta)->required(),
     "Kmer FASTA file");

    umbrella.add(req).add(opt);

    po::store(po::command_line_parser(argc, argv).
              options(umbrella).run(), vm);

    if (vm.count("help") != 0 || argc == 1)
    {
      std::cerr << umbrella << '\n';
      return 1;
    }

    vm.notify();

    uint64_t k = get_k(opts.kmer_fasta);
    std::cerr << "K: " << k << '\n';
    k_map_t kmap = build_reference(opts.genome_fasta, k);
    map(opts.kmer_fasta, kmap);
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << '\n';
    return 1;
  }
}