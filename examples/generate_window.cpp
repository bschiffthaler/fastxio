#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>

int main(int argc, char ** argv)
{

  // Initialize a Record as DNA
  FASTX::Record R("ACGTAGCTAGCTA", "Test", DNA_SEQTYPE);

  // Print 3 bases every 3 bases
  for(auto k : R.window(3,3) )
    {
      std::cout << k << std::endl;
    }

  // Do the same with a window walker
  for(auto w = FASTX::window_walker(3,3,R); ! w.end(); w++)
    {
      std::cout << (*w) << std::endl;
    }

  // Include the final incomplete window
  for(auto w = FASTX::window_walker(3,3,R,true); ! w.end(); w++)
    {
      std::cout << (*w) << std::endl;
    }
  
  return 0;
}
