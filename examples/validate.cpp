#include <iostream>
#include <fastxio_common.h>
#include <fastxio_record.h>
#include <fastxio_reader.h>

int main(int argc, char ** argv)
{


  FASTX::Reader R("100k.fastq", DNA_SEQTYPE);

  std::ios_base::sync_with_stdio(false);
  
  while(R.peek() != EOF)
    {
      FASTX::Record r = R.next();
      if(! r.validate())
	{
	  std::cerr << "Record " << r.get_id() << " could not be validated\n";
	}
    }
    
  return 0;
}
