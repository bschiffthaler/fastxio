#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA with an ambiguous character
  FASTX::Record R("CCGGACTACHVGGGTWTCTAAT", "Primer", DNA_SEQTYPE);

  // Print all possibilities to stdout
  for(auto r : R.enumerate_iupac())
    {
      std::cout << r << std::endl;
    }

  // Initialize a Record as RNA with an ambiguous character
  FASTX::Record R2("NUG", "Test", RNA_SEQTYPE);

  // Print all possibilities to stdout
  for(auto r : R2.enumerate_iupac())
    {
      std::cout << r << std::endl;
    }
  
  return 0;
}
