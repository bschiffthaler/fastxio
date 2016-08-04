#include <iostream>
#include <algorithm> //sort
#include <fastxio_common.h>
#include <fastxio_auxiliary.h> //scan_phred
#include <fastxio_record.h>
#include <fastxio_reader.h>

// A quick and dirty function to calculate percentiles, min, max
std::vector<float> fivenum(std::vector<float> input)
{
  float min = input[0];
  float max = input[input.size() - 1];
  float p25 = input[static_cast<int>(
		      static_cast<float>(
			input.size()) * 0.25 + 1)];
  float p75 = input[static_cast<int>(
		      static_cast<float>(
		        input.size()) * 0.75 + 1)];
  float p50 = input[static_cast<int>(
		      static_cast<float>(
			input.size()) * 0.5 + 1)];
  return std::vector<float>{min, p25, p50, p75, max};
}

int main(int argc, char ** argv)
{

  // Get the PHRED offset
  unsigned short offset = FASTX::scan_phred("p33.fa");
  
  // Open a file as a reader
  FASTX::Reader R("p33.fa", DNA_SEQTYPE);

  std::vector<float> quals; // Aggregate quals
  
  // Loop through records in the file
  while(R.peek() != EOF)
    {
      // Get the next FASTQ record
      FASTX::Record x = R.next();

      unsigned short qsum = 0; // sum of qualities per read
      for(unsigned short q : x.get_numeric_qual() )
	{
	  qsum += (q - offset);
	}
      // Save the mean quality of the read
      quals.push_back( static_cast<float>(qsum) /
		       static_cast<float>(x.size() ));
    }

  // Sort the data
  std::sort(quals.begin(), quals.end());

  // Calculate stats
  std::vector<float> f = fivenum(quals);

  std::cout <<
    "Min: " << f[0] << std::endl <<
    "25%: " << f[1] << std::endl <<
    "50%: " << f[2] << std::endl <<
    "75%: " << f[3] << std::endl <<
    "Max: " << f[4] << std::endl;
  
  return 0;
}
