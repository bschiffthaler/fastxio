#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA
  FASTX::Record R("ACGTAGCTAGCTA", "Test", DNA_SEQTYPE);

  // Print 3 bases every 2 bases
  for(auto k : R.window(3,2) )
    {
      std::cout << k << std::endl;
    }
  
  return 0;
}
