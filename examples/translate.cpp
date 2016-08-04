#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA
  FASTX::Record R("ACGTAGCTAGCTA", "Test", DNA_SEQTYPE);

  // Print all ORFs
  for(auto k : R.translate() )
    {
      std::cout << k << std::endl;
    }

  // Print only the first ORF
  std::cout << R.translate(0) << std::endl;
  
  return 0;
}
