#include <iostream>
#include <fastxio_common.h>
#include <fastxio_reader.h>
#include <fastxio_record.h>
#include <string>

std::string usage = "Usage: "
  "faux_matepair <fastq_file> <insert size> \n\n"
  "Generate a faux mate pair library from a FASTQ file.";

int main(int argc, char ** argv)
{

  if(argc != 3)
    {
      std::cerr << usage << '\n';
      return 1;
    }

  std::string fastq_file = argv[1];
  size_t insert_size = std::stoul(argv[2]);
  size_t fragment_size = insert_size + 200;
  
  // Initialize a Record as DNA
  FASTX::Reader R(fastq_file.c_str(), DNA_SEQTYPE);

  while(R.peek() != EOF)
    {
      FASTX::Record r = R.next();
      for(auto w = FASTX::window_walker(fragment_size, 10, r); ! w.end(); w++)
	{
	  FASTX::Record R1 = (*w).subseq(0,99);
	  FASTX::Record R2 = (*w).subseq(fragment_size - 100, fragment_size);
	  std::cout << R1 << '\n' << R2 << '\n';
	} 
    }
  
  return 0;
}
