#include <iostream>
#include <string>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // First, create a long string for the sequence (900nt)
  std::string seq;
  for(unsigned int i = 0; i < 300; i++)
    {
      seq += "TAA";
    }

  // Initialize a Record as DNA
  FASTX::Record R(seq, "Test", DNA_SEQTYPE);

  // Print wrapped to std::cout
  std::cout << FASTX::Wrap(R) << std::endl;

  // The column width can be specified
  std::cout << FASTX::Wrap(R, 42) << std::endl;
  
  return 0;
}
