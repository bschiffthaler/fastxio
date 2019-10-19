#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <fastxio_gff.h>
#include <fastxio_record.h>
#include <fastxio_reader.h>

struct param
{
  std::string genome_fasta;
  std::string genome_gff;
  uint64_t upstream;
  uint64_t downstream;
  uint64_t minlen;
  bool avoid_clash;
  bool ignore_strand;
};

class GeneBlock
{
public:
  GeneBlock(FASTX::GffRecord const & xgene) :
    found_primary(false), gene(xgene), primary(0) {}

  FASTX::GffRecord const & get_primary();
  FASTX::GffRecord const & get_primary_unsafe() const;
  void find_primary();
  //
  bool operator<(GeneBlock const & other) const {
    return gene < other.gene;
  }
  //
  FASTX::GffRecord gene;
  uint64_t primary;
  bool found_primary;
  mutable std::vector<FASTX::GffRecord> mRNAs;

};

FASTX::GffRecord const & GeneBlock::get_primary()
{
  if (! found_primary)
  {
    find_primary();
  }
  return mRNAs[primary];
}

FASTX::GffRecord const & GeneBlock::get_primary_unsafe() const
{
  return mRNAs[primary];
}


void GeneBlock::find_primary()
{
  uint64_t i = 0;
  uint64_t ir = 0;
  uint64_t best = 0;
  for (auto const & mrna : mRNAs)
  {
    uint64_t start = mrna.start;
    uint64_t end = mrna.end;
    if (start > end)
    {
      std::swap(start, end);
    }
    if (end - start > best)
    {
      best = end - start;
      ir = i;
    }
    i++;
  }
  primary = ir;
  found_primary = true;
}

void get_reg_regions(
  std::vector<GeneBlock> const & blocks,
  FASTX::Gff const & gff,
  std::unordered_map<std::string, FASTX::Record> const & fasta_dict,
  param const & opts,
  bool is_minus)
{
  std::string prev_chrom("");
  for (uint64_t i = 0; i < blocks.size(); i++)
  {
    std::string cur_chrom = blocks[i].gene.seqid;
    FASTX::SequenceRegion seqreg = gff.get_seqregion(cur_chrom);

    uint64_t downstream = opts.downstream;
    uint64_t upstream = opts.upstream;

    if (is_minus)
    {
      std::swap(downstream, upstream);
    }

    try
    {
      fasta_dict.at(cur_chrom);
    }
    catch (std::exception& e)
    {
      throw std::runtime_error("Could not find chromosome key " + cur_chrom +
                               " for gene " +
                               blocks[i].gene.attributes.at("ID") +
                               " in provided FASTX file.");
    }
    FASTX::Record const & fasta_chr = fasta_dict.at(cur_chrom);
    // -----------===========------------
    // ^ left_start         ^ right_start, end
    //            ^ left_end, start
    //                                  ^ right_end
    uint64_t start, end, left_start, left_end, right_start, right_end;

    FASTX::GffRecord const & primary = blocks[i].get_primary_unsafe();
    start = primary.start;
    end = primary.end;
    // Work on + strand
    if (start > end)
    {
      std::swap(start, end);
    }
    left_end = start;
    right_start = end;

    bool no_downstream = false;
    bool no_upstream = false;

    // Check both upstream and downstream for sequence region boundaries
    if (upstream >= start)
    {
      left_start = 1;
    }
    else
    {
      left_start = start - opts.upstream;
    }
    if (end + downstream > fasta_chr.size())
    {
      right_end = fasta_chr.size();
    }
    else
    {
      right_end = end + downstream;
    }

    if (opts.avoid_clash)
    {
      // If we are on a new chromosome, we don't need to check for a previous
      // conflicting gene, but we need to check next
      if (prev_chrom != cur_chrom)
      {
        if ((i + 1) < blocks.size())
        {
          FASTX::GffRecord const & next = blocks[i + 1].get_primary_unsafe();
          uint64_t ns = next.start;
          if (ns > next.end)
          {
            ns = next.end;
          }
          if (next.seqid == cur_chrom)
          {
            if (ns < end)
            { // disable downstream output if we are in another gene
              no_downstream = true;
            }
            if (ns < right_end)
            {
              right_end = ns - 1;
            }
          }
        }
      }
      // Else we check both ends for conflicting genes
      else
      { // previous gene is on the same chromosome
        if ((i + 1) < blocks.size())
        {
          FASTX::GffRecord const & next = blocks[i + 1].get_primary_unsafe();
          uint64_t ns = next.start;
          if (ns > next.end)
          {
            ns = next.end;
          }
          if (next.seqid == cur_chrom)
          { // next gene is on the same chromosome
            if (ns < end)
            { // disable downstream output if we are in another gene
              no_downstream = true;
            }
            if (ns < right_end)
            {
              right_end = ns - 1;
            }
          }
        }
        FASTX::GffRecord const & prev = blocks[i - 1].get_primary_unsafe();
        uint64_t pe = prev.end;
        if (pe < prev.start)
        {
          pe = prev.start;
        }
        if (pe > start)
        { // disable upstream output if we are in another gene
          no_upstream = true;
        }
        if (pe > left_start)
        {
          left_start = pe + 1;
        }
      }
    }

    // Check minimum region length
    if (opts.minlen > 0)
    {
      if (left_end - left_start < opts.minlen)
      {
        no_upstream = true;
      }
      if (right_end - right_start < opts.minlen)
      {
        no_downstream = true;
      }
    }
    if (! no_upstream)
    {
      FASTX::Record sub = fasta_chr.subseq(left_start - 1, left_end - 1);
      sub.set_id(sub.get_id() + 
                 (is_minus ? " type:downstream " : " type:upstream ") +
                 "strand:" + blocks[i].gene.strand + " " +
                 blocks[i].gene.attributes.at("ID"));
      std::cout
          << FASTX::Wrap(sub)
          << '\n';
    }
    if (! no_downstream)
    {
      FASTX::Record sub = fasta_chr.subseq(right_start - 1, right_end - 1);
      sub.set_id(sub.get_id() +
                 (is_minus ? " type:upstream " : " type:downstream ") +
                 "strand:" + blocks[i].gene.strand + " " +
                 blocks[i].gene.attributes.at("ID"));
      std::cout
          << FASTX::Wrap(sub)
          << '\n';
    }
    prev_chrom = cur_chrom;
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
    ("help,h", "Show this help message")
    ("upstream,u", po::value<uint64_t>(&opts.upstream)->default_value(1000),
     "How far to extract upstream")
    ("downstream,d", po::value<uint64_t>(&opts.downstream)->default_value(1000),
     "How far to search downstream")
    ("minlen,m", po::value<uint64_t>(&opts.minlen)->default_value(50),
     "Minimum length of promoter region to be output")
    ("avoid_clash,a", po::bool_switch(&opts.avoid_clash),
     "Truncate regions if they clash with the beginning or end of another gene")
    ("ignore_strand,i", po::bool_switch(&opts.ignore_strand),
     "Ignore strand information.");

    po::options_description req("Required");
    req.add_options()
    ("fasta,f", po::value<std::string>(&opts.genome_fasta)->required(),
     "Genome FASTA file")
    ("gff,g", po::value<std::string>(&opts.genome_gff)->required(),
     "Genome GFF file");

    umbrella.add(req).add(opt);

    po::store(po::command_line_parser(argc, argv).
              options(umbrella).run(), vm);

    if (vm.count("help") != 0 || argc == 1)
    {
      std::cerr << umbrella << '\n';
      return 1;
    }

    vm.notify();

    std::unordered_map<std::string, FASTX::Record> fasta_dict;
    FASTX::Reader Reader(opts.genome_fasta.c_str(), DNA_SEQTYPE);
    while (Reader.peek() != EOF)
    {
      FASTX::Record rec = Reader.next();
      fasta_dict[rec.get_canonical_id()] = rec;
    }

    FASTX::Gff gff(opts.genome_gff);

    std::unordered_map<std::string, GeneBlock> gff_genes;
    // First get all genes
    for (uint64_t i = 0; i < gff.records().size(); i++)
    {
      FASTX::GffRecord const & current = gff.records()[i];
      if (current.type == "gene")
      {
        std::string const key(current.attributes.at("ID"));
        gff_genes.emplace(key, GeneBlock(current));
      }
    }
    // Now get all mRNAs
    for (uint64_t i = 0; i < gff.records().size(); i++)
    {
      FASTX::GffRecord const & current = gff.records()[i];
      if (current.type == "mRNA")
      {
        std::string const parent_id(current.attributes.at("Parent"));
        gff_genes.at(parent_id).mRNAs.push_back(current);
      }
    }
    // Identify primary mRNA
    for (auto & gb : gff_genes)
    {
      gb.second.find_primary();
    }
    // Get blocks into a vector
    std::vector<GeneBlock> all_blocks;
    all_blocks.reserve(gff_genes.size());
    for (auto & key_value : gff_genes)
    {
      all_blocks.push_back(key_value.second);
    }
    std::sort(all_blocks.begin(), all_blocks.end());

    if (opts.ignore_strand)
    {
      get_reg_regions(all_blocks, gff, fasta_dict, opts, false);
    }
    else
    {
      std::vector<GeneBlock> plus;
      std::vector<GeneBlock> minus;
      for (auto & g : all_blocks)
      {
        if (g.gene.strand == '+')
        {
          plus.push_back(g);
        }
        else if (g.gene.strand == '-')
        {
          minus.push_back(g);
        }
        else
        {
          std::cerr << "[WARNING]: Gene " << g.gene.attributes.at("ID") 
          << " is not marked to be on the + or - strand and will be ignored.\n"; 
        }
      }
      get_reg_regions(plus, gff, fasta_dict, opts, false);
      get_reg_regions(minus, gff, fasta_dict, opts, true);
    }

  }
  catch (std::exception& e)
  {
    std::cerr << "[ERROR]: " << e.what() << '\n';
    return 1;
  }
  return 0;
}