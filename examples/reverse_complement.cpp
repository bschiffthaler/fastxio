#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA
  FASTX::Record R("ACGTAGCTAGCTA", "Test", DNA_SEQTYPE);

  // Print the original
  std::cout << R << std::endl;
  
  // Print the reverse complement
  std::cout << R.rc() << std::endl;
  
  return 0;
}
