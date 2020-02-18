// Copyright (2018) <BAstian Schiffthaler>
#include <omp.h>
#include <fastxio_reader.h>
#include <fastxio_record.h>
#include <fastxio_minhash.h>
#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, char ** argv)
{
  try
  {
    uint16_t threads;  // number of threads
    uint32_t nhash;  // number of hashes
    uint32_t k;  // length of k
    std::string ref;  // reference sequence file
    std::string target;  // target sequence file
    uint32_t batch_size;  // vector buffer size for parallel processing
    double sim_cutoff;  // Similarity cutoff to output filtered data
    std::string cont_out;
    std::string clean_out;
    bool print_stats;

    po::options_description umbrella;
    po::options_description opt("Options");
    opt.add_options()
    ("batch-size,b", po::value<uint32_t>(&batch_size)->default_value(1024),
     "Number of records in buffer for parallel processing")
    ("contaminant,c", po::value<std::string>(&cont_out)->default_value(""),
     "Output file for contaminant sequences")
    ("help,h", "Show this help message")
    ("k-length,k", po::value<uint32_t>(&k)->default_value(13),
     "K-mer length")
    ("min-similarity,m", po::value<double>(&sim_cutoff)->default_value(-1),
     "Minimal similarity for a read to be considered a match")
    ("nhash,n", po::value<uint32_t>(&nhash)->default_value(100),
     "Number of min hashes to keep")
    ("out,o", po::value<std::string>(&clean_out)->default_value(""),
     "Output file for clean sequences")
    ("print-stats,s", po::bool_switch(&print_stats)->default_value(false))
    ("threads,t", po::value<uint16_t>(&threads)->default_value(1),
     "Number of parallel threads");

    po::options_description req("Required");
    req.add_options()
    ("ref", po::value<std::string>(&ref)->required(),
     "Input reference sequence")
    ("target", po::value<std::string>(&target)->required(),
     "Input target sequence");

    umbrella.add(opt).add(req);

    po::positional_options_description p;
    p.add("ref", 1);
    p.add("target", 1);

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
    FASTX::MinHash hash(nhash, k);
    while (ref_reader.peek() != EOF)
    {
      FASTX::Record r = ref_reader.next();
      hash.add(r);
    }

    // Prepare to read target sequences into a vector
    // as batches
    FASTX::Reader test(target.c_str(), DNA_SEQTYPE);
    // std::vector<FASTX::Record> tvec;
    // tvec.reserve(batch_size);

    // Prepare output files
    std::ofstream clean_ost;
    std::ofstream cont_ost;

    if (clean_out != "")
      clean_ost.open(clean_out.c_str(), std::ios::out);
    if (cont_out != "")
      cont_ost.open(cont_out.c_str(), std::ios::out);

    omp_lock_t cplock;
    omp_init_lock(&cplock);
    omp_set_num_threads(threads);

    #pragma omp parallel
    {
      #pragma omp single
      {
        // Read all target sequences
        while (test.peek() != EOF)
        {
          // Add to batch
          // tvec.push_back(test.next());
          // if (tvec.size() == batch_size || test.peek() == EOF)
          // {
          //   // Process batches in parallel
          //   #pragma omp parallel for
          //   for (uint64_t i = 0; i < tvec.size(); i++)
          //   {
          auto seq = test.next();
          #pragma omp task default(shared) firstprivate(seq)
          {
            // Get reference record with best similary to target record
            auto sim = hash.max_similarity(seq);
            // Write output only one thread at a time
            // target_id  reference_id  hits  nhash_a nhash_b jaccard_similarity
            omp_set_lock(&cplock);
            if (print_stats)
            {
              // Write match statistics
              std::cout
                  << seq.get_id() << '\t'
              << hash.id(sim.idx) << '\t'
              << sim.hits << '\t'
              << sim.asize << '\t'
              << sim.bsize << '\t'
              << sim.ji << '\n';
            }
            // Write FASTX
            if (sim.ji > sim_cutoff)
            {
              if (cont_out != "")
              {
                cont_ost
                    << seq.get_id() << " || " << hash.id(sim.idx) << '\n'
                    << seq.get_seq() << '\n';
              }
            }
            else if (clean_out != "")
            {
              clean_ost << seq << '\n';
            }

            omp_unset_lock(&cplock);
          }
          // Reset batch vector
          //   tvec.clear();
          //   tvec.shrink_to_fit();
          // }
        }
        #pragma omp taskwait
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
