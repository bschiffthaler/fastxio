#include <fastxio_minhash.h>
#include <fastxio_record.h>
#include <MurmurHash3.h>
#include <cstdint>
#include <vector>
#include <unordered_set>

#define SEED 314159265

namespace FASTX {

extern GData global;

uint32_t set_isec(std::unordered_set<uint32_t>& lhs,
                  std::unordered_set<uint32_t>& rhs)
{
  uint32_t hits = 0;
  for (auto& item : lhs)
  {
    if (rhs.find(item) != rhs.end())
      hits++;
  }
  return hits;
}

MinHash::MinHash(uint32_t j, uint32_t k) : _j(j), _k(k)
{}

void MinHash::add(const Record& rec, bool complement)
{
  std::vector<uint32_t> hashes;
  hashes.reserve(rec.size() - _k + 1);
  int _ki = _k;
  for (uint64_t i = 0; i < (rec.size() - _k); i++)
  {
    uint32_t hash[4]; // We get 128 bits from MurmurHash3_x64_128
    MurmurHash3_x64_128(&rec._seq[i], _ki, SEED, &hash);
    hashes.push_back(hash[0]); // Keep only 32 bits
  }
  std::sort(hashes.begin(), hashes.end());
  std::unordered_set<uint32_t> hashset;
  for (auto& h : hashes)
  {
    hashset.insert(h);
    if (hashset.size() == _j) break;
  }
  _hashes.push_back(hashset);
  _ids.push_back(rec._id);
  if (complement)
  {
    Record tmp = !rec;   
    this->add(tmp, false);
  }
}

minhash_sim_t MinHash::max_similarity(const Record& rec)
{
  std::vector<uint32_t> hashes;
  hashes.reserve(rec.size() - _k + 1);
  for (uint64_t i = 0; i < (rec.size() - _k); ++i)
  {
    uint32_t hash[4]; // We get 128 bits from MurmurHash3_x64_128
    MurmurHash3_x64_128(&rec._seq[i], _k, SEED, &hash);
    hashes.push_back(hash[0]); // Keep only 32 bits
  }
  std::sort(hashes.begin(), hashes.end());
  std::unordered_set<uint32_t> hashset;
  for (auto& h : hashes)
  {
    hashset.insert(h);
    if (hashset.size() == _j) break;
  }

  minhash_sim_t res;

  res.ji = -std::numeric_limits<double>::infinity();
  uint32_t ctr = 0;
  for(auto& comp : _hashes)
  {
    uint32_t hits = set_isec(hashset, comp);
    double hitsd = hits;
    double asize = hashset.size();
    double bsize = comp.size();
    double sim = hitsd / (asize + bsize - hitsd);
    if (sim > res.ji)
    {
      res.asize = hashset.size();
      res.hits = hits;
      res.bsize = comp.size();
      res.ji = sim;
      res.idx = ctr;
    }
    ctr++;
  }
  return res;
}

}