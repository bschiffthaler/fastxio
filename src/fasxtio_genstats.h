#ifndef _FASTX_IO_GENOME_STATS_H_
#define _FASTX_IO_GENOME_STATS_H_

#include <string>
#include <memory>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <matrix.h>
#include <fastxio_common.h>
#include <fastxio_nuc_frequency.h>

namespace FASTX {

/**
 * @brief Compute and store genome statistics.
 *
 */
class GenomeStats {
public:


private:
  NucFrequency _nuc_freq;
  size_t _length_total;
  size_t _n50;
  size_t _l50;
  size_t _shortest;
  size_t _longest;
};

}
#endif
