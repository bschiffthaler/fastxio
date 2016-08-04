#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#ifndef NO_ERROR_CHECKING
  #include <cerrno>
  #include <stdexcept>
#endif

#include <gzstream.h>
#include <bz2stream.h>
#include <matrix.h>
#include <set>
#include <map>
#include <fastxio_record.h>
#include <fastxio_common.h>
#include <fastxio_nuc_frequency.h>
#include <fastxio_auxiliary.h>

namespace FASTX {

  extern const std::set<char> G_nuc_alphabet;
  extern const std::set<char> G_aa_alphabet;
  extern const std::map<char, char> G_rc;
  extern const std::map<char, std::vector<char> > G_enum_iupac_dna;
  extern const std::map<char, std::vector<char> > G_enum_iupac_rna;
  extern const std::map<std::string, char> G_codon_to_protein;
  
  //Istream constructor
  Record::Record(std::istream& input, char seqtype = DNA_SEQTYPE)
    {
      std::getline(input, _id);
      if(_id.front() == '>')
	_type = (FASTA_TYPE | seqtype);
      else if(_id.front() == '@')
	_type = (FASTQ_TYPE | seqtype);
      else
	throw std::runtime_error("Could not determine format: " + _id);
      _id.erase(0,1);
      std::getline(input, _seq);
      if(_type & FASTA_TYPE)
	{
	  while(is_sequence_char(input.peek(), seqtype))
	    {
	      std::string app;
	      std::getline(input, app);
	      _seq.append(app);
	    }
	}
      if(_type & FASTQ_TYPE)
	{
	  input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	  std::getline(input, _qual); 
	}
      #ifndef NO_ERROR_CHECKING
      this->validate();
      #endif
    }


  //FASTA sequence constructor
  Record::Record(const std::string& seq, const std::string& id,
		 char seqtype = DNA_SEQTYPE) :
    _seq(seq), _id(id), _type(FASTA_TYPE | seqtype)
  {
#ifndef NO_ERROR_CHECKING
    this->validate();
#endif
  }

  //FASTQ sequence constructor
  Record::Record(const std::string& seq, const std::string& id,
		 const std::string& qual, char seqtype = DNA_SEQTYPE) :
    _seq(seq), _id(id), _qual(qual), _type(FASTQ_TYPE | seqtype)
  {
#ifndef NO_ERROR_CHECKING
    this->validate();
#endif
  }

  bool Record::validate(void) const
  {
    bool ret = true;
    // FASTQ
    if(_type & FASTQ_TYPE)
      {
	if(_qual.size() != _seq.size())
	  {
	    throw std::runtime_error("Qual and sequence are not the same length for: " + _id);
	    ret = false;
	  }
      }
    if((_type & DNA_SEQTYPE) || (_type & RNA_SEQTYPE))
      {
	for(char c : _seq)
	  {
	    if(G_nuc_alphabet.find(c) == G_nuc_alphabet.end())
	      {
		throw std::runtime_error("Unknown character " + std::to_string(c) + " in sequence: " + _id);
		ret = false;
	      }
	  }
      }
    if(_type & AA_SEQTYPE)
      {
	for(char c : _seq)
	  {
	    if(G_aa_alphabet.find(c) == G_aa_alphabet.end())
	      {
		throw std::runtime_error("Unknown character " + std::to_string(c) + " in sequence: " + _id);
		ret = false;
	      }
	  }
      }
    if(_type & FASTQ_TYPE)
      {
	for(char c : _qual)
	  {
	    if(c < 33 || c > 104)
	      {
		throw std::runtime_error("Impossible quality score " + std::to_string(c) + " in sequence: " + _id);
		ret = false;
	      }
	  }
      }
    return ret;
  }
  
  //Print output
  std::ostream& operator<<(std::ostream& outstream, const Record& rec)
  {
    if(rec._type & FASTQ_TYPE)
      outstream << '@' << rec._id << '\n' << rec._seq << "\n+\n" << rec._qual;
    if(rec._type & FASTA_TYPE)
      outstream << '>' << rec._id << '\n' << rec._seq;
    return outstream;
  }

  
  std::vector<unsigned short> Record::get_numeric_qual(void) const
  {
#ifndef NO_ERROR_CHECKING
    if(_type & FASTA_TYPE)
      throw std::runtime_error("FASTA files do not have quality values");
#endif
    std::vector<unsigned short> ret;
    for(unsigned char c : _qual) 
      ret.push_back(c);
    return(ret);
  }
  
  // Translate DNA/RNA to protein
  Record Record::translate(unsigned short frame)
  {
#ifndef NO_ERROR_CHECKING
    if(_type & AA_SEQTYPE)
      throw std::runtime_error("Cannot translate amino acid sequence");
#endif
    std::stringstream ss(_seq);
    std::string res;
    ss.ignore(frame, EOF);
    unsigned short count = 0;
    std::string tmp;
    while(ss.good())
      {
        if(count == 3)
          {
            count = 0;
            auto target = G_codon_to_protein.find(tmp);
            if(target == G_codon_to_protein.end())
              {
                res += 'X';
              }
            else
              {
                res += target->second; 
              }
            tmp.clear();
          }
        tmp += ss.get();
        count++;
      }
    return Record(res, _id + " ORF" + std::to_string(frame), AA_SEQTYPE);
  }

  std::vector<Record>  Record::translate(void)
  {
    std::vector<Record> res;
    res.push_back( this->translate(0) );
    res.push_back( this->translate(1) );
    res.push_back( this->translate(2) );
    return res;
  }

  // Extract a subsequence
  Record Record::subseq(length_t start, length_t stop) const
  {
    std::string seq = _seq.substr(start, (stop - start) + 1);
    std::string qual;
    if(_type & FASTQ_TYPE)
      qual = _qual.substr(start, (stop - start) + 1);
    std::string id = _id + " " + std::to_string(start) + "-" + 
      std::to_string(stop);

    if(_type & FASTA_TYPE)
      return Record(seq, id);
    else
      return Record(seq, id, qual);
  }

  // Create kmers along a sequence
  std::vector<Record> Record::kmer(length_t k) const
  {
    length_t start = 0;
    length_t stop = k - 1;
    std::vector<Record> res;
    while(stop < _seq.length())
      {
	res.push_back(this->subseq(start++, stop++));
      }
    return res;
  }

  // Create a sliding window along a sequence
  std::vector<Record> Record::window(length_t width, length_t increment) const
  {
    length_t start = 0;
    length_t stop = width - 1;
    std::vector<Record> res;
    while(stop < _seq.length())
      {
	res.push_back(this->subseq(start, stop));
	start += increment;
	stop += increment;
	if(stop >= _seq.length())
	  res.push_back(this->subseq(start, _seq.length() - 1));
      }
    return res;
  }

  // Get nucleotide frequency of a record
  NucFrequency Record::count_freq(void)
  {
    NucFrequency ret;
    ret.add(_seq);
    return ret;
  }

  // Reverse complement
  Record Record::rc() const
  {
    std::string seq(_seq.rbegin(), _seq.rend());
    for(auto it = seq.begin(); it != seq.end(); it++)
      {
	*it = G_rc.at(*it);
      }
    if(_type & FASTQ_TYPE)
      {
        std::string qual(_qual.rbegin(), _qual.rend());
        return Record(seq, _id + " RC", qual);
      }
    else
      {
        return Record(seq, _id + " RC");
      }
  }

  std::set<Record> Record::enumerate_iupac(void)
  {
    std::set<Record> res;
    if(this->_type & DNA_SEQTYPE)
      {
	std::set<char> unambiguous_nuc = {'A', 'C', 'G', 'T'};
	recursive_iupac_enum(res, *this, G_enum_iupac_dna, unambiguous_nuc);
      }
    else if(this->_type & RNA_SEQTYPE)
      {
	std::set<char> unambiguous_nuc = {'A', 'C', 'G', 'U'};
	recursive_iupac_enum(res, *this, G_enum_iupac_rna, unambiguous_nuc);
      }
    return res;
  }

  Wrap::Wrap(const Record& rec)
  {
#ifndef NO_ERROR_CHECKING
    if(rec._type & FASTQ_TYPE)
      throw std::runtime_error("Cannot line-wrap a FASTQ record.");
#endif
    _id = std::make_shared<const std::string>(rec._id);
    _seq = std::make_shared<const std::string>(rec._seq);
  }

  Wrap::Wrap(const Record& rec, unsigned int width)
  {
#ifndef NO_ERROR_CHECKING
    if(rec._type & FASTQ_TYPE)
      throw std::runtime_error("Cannot line-wrap a FASTQ record.");
#endif
    _id = std::make_shared<const std::string>(rec._id);
    _seq = std::make_shared<const std::string>(rec._seq);
    _width = width;
  }

  //Print wrapping at N chars
  std::ostream& Wrap::operator()(std::ostream& outstream) const
  {
    outstream << '>' << *_id << '\n';
    unsigned short count = 0;
    size_t len = _seq->length();
    size_t pos = 0;
    for(char c : *_seq)
      {
	if(count == _width && pos != (len-1) )
	  {
	    outstream << '\n';
	    count = 0;
	  }
	outstream << c;
	count++;
	pos++;
      }
    return outstream;
  }

  std::ostream& operator<<(std::ostream& outstream, Wrap rec)
  {
    return rec(outstream);
  }
  
};
