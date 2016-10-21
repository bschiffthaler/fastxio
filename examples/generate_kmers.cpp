#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA
  FASTX::Record R("ACGTAGCTAGCTA", "Test", DNA_SEQTYPE);

  // Generate a vector containing kmers of length 5 and print them
  for(auto k : R.kmer(5) )
    {
      std::cout << k << std::endl;
    }

  // Do the same using a kmer walker
  for(auto k = FASTX::kmer_walker(5, R); ! k.end(); k++)
    {
      std::cout << (*k) << std::endl;
    }
  
  return 0;
}
