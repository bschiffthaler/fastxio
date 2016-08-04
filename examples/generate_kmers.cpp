#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA
  FASTX::Record R("ACGTAGCTAGCTA", "Test", DNA_SEQTYPE);

  // Print all k-mers of length 3
  for(auto k : R.kmer(3) )
    {
      std::cout << k << std::endl;
    }
  
  return 0;
}
