# FASTX-IO

_No BS FASTA/Q input/output and other operations_

This library is meant as a simple, easy-to-use and effective means of reading, writing and manipulating [FASTA] and [FASTQ] files. The emphasis of the library lies on a balance of speed and correctness, therefore it aims to be IUPAC compliant, supporting ambiguous nucleotide letters (e.g. `W`) and DNA, RNA, as well as amino acid containing data.

It is meant as building blocks for larger applications.

## Features
_See also the examples at the end of this document_

* Easy I/O from files
* Automatic detecting and reading of `gzip` and `bzip2` compressed files
* Built-in support for many common operations
    * Simple generation of sub-records:
	    * k-mers
	    * sliding windows
	    * sub-sequences
    * Translating DNA/RNA to amino acids
    * Reverse complementing
    * Tabulating counts of nucleotide frequencies
    * Fast enumeration of all possible sequences from ambiguous ones (e.g. primers)
* Automatic validation of records
* Optional compilation without record validation for extra speed


## Speed

`fastxio` performs minimal validation of the records and makes use of many speed efficient data structures to perform operations. Record validation can be turned off completely during compilation, though only advisable if the data has been validated already before (e.g. by `fastQValidator`).

## Building the library

`cmake` is required to build the library, `doxygen` to build the documentation.

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

In order to build the library without record validation, pass the `NO_CHECKING` flag to `cmake`

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DNO_CHECKING=ON ..
make
```


## Examples
### Counting the frequencies of all records and calculating GC%

Calculating the GC% of a record can often be useful to guess at the locality of a sequency (e.g. rDNA) or at the origin (fungal contamination). `fastxio` employs a nucleotide counting class:

```
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
``` 

**Output:**

```
A	10024
C	4901
G	5106
N	110
T	9959
0.33%
```

### Enumerating all possible real sequences from an ambiguous primer

Often, when designing primers, these can contain ambiguous nucleotides, but many detection and trimming programs do not support anything but A,C,G and T. Therefore, it can be useful to have a way of generating all (non-redundant) possibilities.

```
// Initialize a Record as DNA with an ambiguous character
FASTX::Record R("CCGGACTACHVGGGTWTCTAAT", "Primer", DNA_SEQTYPE);

// Print all possibilities to stdout
for(auto r : R.enumerate_iupac())
  {
    std::cout << r << std::endl;
  }
```

**Output:**

```
>Primer_10A_11A_16A
CCGGACTACAAGGGTATCTAAT
>Primer_10A_11A_16T
CCGGACTACAAGGGTTTCTAAT
>Primer_10A_11C_16A
CCGGACTACACGGGTATCTAAT
>Primer_10A_11C_16T
CCGGACTACACGGGTTTCTAAT
>Primer_10A_11G_16A
CCGGACTACAGGGGTATCTAAT
>Primer_10A_11G_16T
CCGGACTACAGGGGTTTCTAAT
>Primer_10C_11A_16A
CCGGACTACCAGGGTATCTAAT
>Primer_10C_11A_16T
CCGGACTACCAGGGTTTCTAAT
>Primer_10C_11C_16A
CCGGACTACCCGGGTATCTAAT
>Primer_10C_11C_16T
CCGGACTACCCGGGTTTCTAAT
>Primer_10C_11G_16A
CCGGACTACCGGGGTATCTAAT
>Primer_10C_11G_16T
CCGGACTACCGGGGTTTCTAAT
>Primer_10T_11A_16A
CCGGACTACTAGGGTATCTAAT
>Primer_10T_11A_16T
CCGGACTACTAGGGTTTCTAAT
>Primer_10T_11C_16A
CCGGACTACTCGGGTATCTAAT
>Primer_10T_11C_16T
CCGGACTACTCGGGTTTCTAAT
>Primer_10T_11G_16A
CCGGACTACTGGGGTATCTAAT
>Primer_10T_11G_16T
CCGGACTACTGGGGTTTCTAAT
```

