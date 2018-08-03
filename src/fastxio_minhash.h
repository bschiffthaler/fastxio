#ifndef _FASTX_IO_MINHASH_H_
#define _FASTX_IO_MINHASH_H_

#include <fastxio_record.h>
#include <cstdint>
#include <unordered_set>

namespace FASTX {

struct minhash_sim_t
{
  double ji = 0;
  uint32_t hits = 0;
  uint32_t idx = 0;
  uint32_t asize = 0;
  uint32_t bsize = 0;
};

/**
   * @brief Compute the number of intersects of two unordered sets.
   * 
   * @param lhs First set
   * @param rhs Second set
   * @return The number of items shared between the two sets.
   */
uint32_t set_isec(std::unordered_set<uint32_t>& lhs,
                  std::unordered_set<uint32_t>& rhs);

/** 
   * @brief kmer overlap based on MinHash.
   *
   * This class will produce kmer overlaps using a locality sensitive
   * hashing approach (MinHash). Specifically, the implementation is
   * the single hash function variant.
   */
class MinHash {
public:
  /**
     * @brief Empty constructor. The number of hashes and kmer size must be
     * chosen at construct time.
     *
     * @param j The number of hashes
     * @param k The sequence kmer size
     *
     * @return A MinHash object with the given parameters.
     */
  MinHash(uint32_t j, uint32_t k);
  /**
     * @brief Add a `Record` object to the kmer index.
     *
     * @param r A `Record` object containing the sequence that is to be added to
     * the index.
     * @param complement `true` if the reverse complement of the sequence should
     * also be added to the index.
     *
     */
  void add(const Record& r, bool complement = true);
  /**
     * @brief Get the maximum similarity of a new `Record` object to the index.
     * (best hit)
     *
     * @param r A new `Record` object.
     *
     * @return A pair of values. The first one is the match similarity, the
     * second points to the index of id of the sequence that was the target.
     * Meant to be used with `MinHash::id()`
     */
  minhash_sim_t max_similarity(const Record& r);
  /**
     * @brief Get the ID of an index sequence.
     *
     * @param idx The index that was returned by e.g.: `MinHash::max_similarity()`.
     * 
     * @return The ID of the sequence that matches the index.
     */
  std::string& id(uint32_t idx){return _ids[idx];}
private:
  // Number of hashes
  const uint32_t _j;
  // Kmer length
  const uint32_t _k;
  std::vector<std::unordered_set<uint32_t>> _hashes;
  std::vector<std::string> _ids;
};
}

#endif