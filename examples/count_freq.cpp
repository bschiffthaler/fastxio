#include <iostream>
#include <iomanip> //setprecision
#include <fastxio_common.h>
#include <fastxio_record.h>
#include <fastxio_reader.h>
#include <fastxio_nuc_frequency.h>

int main(int argc, char ** argv)
{

  // Open a file as a reader
  FASTX::Reader R("p33.fa", DNA_SEQTYPE);

  // Create a NucFrequency object to keep track
  FASTX::NucFrequency NF;
  
  // Loop through records in the file
  while(R.peek() != EOF)
    {
      // Get the next FASTQ record
      FASTX::Record x = R.next();
      // Add it to the count
      NF.add(x);
    }

  // Print the result
  std::cout << NF;

  // Calculate GC%
  std::cout << std::setprecision(2) << static_cast<float>( NF['G'] + NF['C'] ) /
    static_cast<float>( NF['G'] + NF['C'] + NF['A'] + NF['T'] ) << "%" << std::endl;

  
  return 0;
}
