#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <cctype>

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
 
  extern GData global; /**< Global variable of translation tables */

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
	  while(isspace(input.peek()))
	    input.ignore(1, '\n');
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
    _seq(seq), _id(id),  _type(FASTA_TYPE | seqtype)
  {
#ifndef NO_ERROR_CHECKING
    this->validate();
#endif
  }

  //FASTQ sequence constructor
  Record::Record(const std::string& seq, const std::string& id,
		 const std::string& qual, char seqtype = DNA_SEQTYPE) :
    _seq(seq), _id(id), _qual(qual),  _type(FASTQ_TYPE | seqtype)
  {
#ifndef NO_ERROR_CHECKING
    this->validate();
#endif
  }

  Record::Record(void) :
    _seq(""), _id(""), _qual(""), _type(NULL_SEQTYPE)
  {
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
	    if(global.nuc_alphabet.find(c) == global.nuc_alphabet.end())
	      {
          std::string errmsg = "Unknown character ";
          errmsg += c;
		throw std::runtime_error(errmsg + " in sequence: " + _id);
		ret = false;
	      }
	  }
      }
    if(_type & AA_SEQTYPE)
      {
	for(char c : _seq)
	  {
	    if(global.aa_alphabet.find(c) == global.aa_alphabet.end())
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
    auto& translation_table = _type & DNA_SEQTYPE ? global.codon_to_protein_dna :
                              global.codon_to_protein_rna;
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
            auto target = translation_table.find(tmp);
            if(target == translation_table.end())
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
  std::vector<Record> Record::window(length_t width, length_t increment,
				     bool include_final) const
  {
    length_t start = 0;
    length_t stop = width - 1;
    std::vector<Record> res;
    while(stop < _seq.length())
      {
	res.push_back(this->subseq(start, stop));
	start += increment;
	stop += increment;
	if(stop >= _seq.length() && include_final)
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
	*it = global.rc.at(*it);
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
	recursive_iupac_enum(res, *this, global.enum_iupac_dna, unambiguous_nuc);
      }
    else if(this->_type & RNA_SEQTYPE)
      {
	std::set<char> unambiguous_nuc = {'A', 'C', 'G', 'U'};
	recursive_iupac_enum(res, *this, global.enum_iupac_rna, unambiguous_nuc);
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

  //Initialize kmer walker at position 0
  kmer_walker::kmer_walker(size_t k, const Record& r) :
    _k(k), _parent(r), _current_pos(0)
  {
    if(_parent.size() >= _k)
      {
	_rec = _parent.subseq(0, _k - 1); // Stop of subseq is inclusive
	_end = false;
      }
    else
      {
	_end = true;
      }
    _begin = false;
  }

  // Initialize a kmer walker at a given position
  kmer_walker::kmer_walker(size_t k, size_t pos, const Record& r) :
    _k(k), _parent(r), _current_pos(pos)
  {
    if(_parent.size() >= (_k + _current_pos))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_end = false;
      }
    else
      {
	_end = true;
      }
    _begin = _current_pos == 0 ? true : false;
  }

  // Move kmer by +1 in pre increment
  kmer_walker& kmer_walker::operator++()
  {
    _current_pos++;
    _begin = false; //by definition
    
    if(_parent.size() >= (_current_pos + _k))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_end = false;
      }
    else
      {
	_end = true;
      }
    
    return *this;
  }

  // Move kmer by +1 in post increment
  kmer_walker kmer_walker::operator++(int)
  {
    kmer_walker ret(*this);
    _current_pos++;
    _begin = false;
    
    if(_parent.size() >= (_current_pos + _k))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_end = false;
      }
    else
      {
	_end = true;
      }
    
    return ret;
  }

  // Move kmer by -1 in pre decrement
  kmer_walker& kmer_walker::operator--()
  {
    _end = _current_pos <= (_parent.size() - 1) ? false : true;
    if(_current_pos == 0)
      {
	_begin = true;
      }
    else
      {
	_current_pos--;
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_begin = false;
      }
    return *this;
  }

  // Move kmer by -1 in post decrement
  kmer_walker kmer_walker::operator--(int)
  {
    _end = _current_pos <= (_parent.size() - 1) ? false : true;
    kmer_walker ret(*this);
    if(_current_pos == 0)
      {
	_begin = true;
      }
    else
      {
	_current_pos--;
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_begin = false;
      }
    return ret;
  }

  // Return reference to kmer
  Record& kmer_walker::operator*()
  {
    return _rec;
  }

  // Return a copy of the kmer
  Record kmer_walker::kmer()
  {
    return _rec;
  }

  // Skip a number of positions
  void kmer_walker::skip(size_t n)
  {
    _current_pos += n;
    if(_parent.size() >= (_current_pos + _k))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_end = false;
      }
    else
      {
	_end = true;
      }
  }

  // Rewind a number of bases
  void kmer_walker::rewind(size_t n)
  {
    if(_current_pos < n)
      {
	_begin = true;
      }
    else
      {
	_current_pos -= n;
	_rec = _parent.subseq(_current_pos, _current_pos + _k - 1);
	_begin = false;
      }
  }

  //Initialize a window walker at position 0
  window_walker::window_walker(size_t ws, size_t increment,
			       const Record& r, bool include_final) :
    _ws(ws), _parent(r), _current_pos(0),
    _increment(increment), _include_final(include_final)
  {
    if(_parent.size() >= _ws)
      {
	_rec = _parent.subseq(0, _ws - 1); // Stop of subseq is inclusive
	_end = false;
      }
    else if(_include_final)
      {
	_rec = _parent.subseq(0, _parent.size() - 1);
	_end = false; _include_final = false;
      }
    else
      {
	_end = true;
      }
    _begin = true;
  }

  // Initialize a window walker at a given position
  window_walker::window_walker(size_t ws, size_t increment,
			       size_t pos, const Record& r,
			       bool include_final) :
    _ws(ws), _parent(r),
    _current_pos(pos), _increment(increment)
  {
    if(_parent.size() <= (_ws + _current_pos))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_end = false;
      }
    else if(_include_final && _current_pos < _parent.size())
      {
	_rec = _parent.subseq(_current_pos, _parent.size() - 1);
	_end = false; _include_final = false;
      }
    else
      {
	_end = true;
      }
    _begin = _current_pos == 0 ? true : false;
  }

  // Move window by +increment in pre increment
  window_walker& window_walker::operator++()
  {
    _current_pos += _increment;
    _begin = false;
    if(_parent.size() >= (_current_pos + _ws))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_end = false;
      }
    else if(_include_final && _current_pos < _parent.size())
      {
	_rec = _parent.subseq(_current_pos, _parent.size() - 1);
	_end = false; _include_final = false;
      }
    else
      {
	_end = true;
      }
    return *this;
  }

  // Move window by +increment in post increment
  window_walker window_walker::operator++(int)
  {
    window_walker ret(*this);
    _current_pos += _increment;
    _begin = false;
    if(_parent.size() >= (_current_pos + _ws))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_end = false;
      }
    else if(_include_final && _current_pos < _parent.size())
      {
	_rec = _parent.subseq(_current_pos, _parent.size() - 1);
	_end = false; _include_final = false;
      }

    else
      {
	_end = true;
      }
    return ret;
  }

  // Move window by -increment in pre decrement
  window_walker& window_walker::operator--()
  {
    if(_current_pos < _increment)
      {
	if(_include_final)
	  {
	    _rec = _parent.subseq(0, _current_pos);
	    _begin = false; _include_final = false;
	  }
	else
	  {
	    _begin = true;
	  }
	_end = _current_pos <= (_parent.size() - 1) ? false : true;
      }
    else
      {
	_current_pos -= _increment;
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_begin = false;
	_end = _current_pos <= (_parent.size() - 1) ? false : true;
      }
    return *this;
  }

  // Move window by -increment in post decrement
  window_walker window_walker::operator--(int)
  {
    window_walker ret(*this);
    if(_current_pos < _increment)
      {
	if(_include_final)
	  {
	    _rec = _parent.subseq(0, _current_pos);
	    _begin = false; _include_final = false;
	  }
	else
	  {
	    _begin = true;
	  }
	_end = _current_pos <= (_parent.size() - 1) ? false : true;
      }
    else
      {
	_current_pos -= _increment;
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_begin = false;
	_end = _current_pos <= (_parent.size() - 1) ? false : true;
      }
    return ret;
  }

  // Return reference to window
  Record& window_walker::operator*()
  {
    return _rec;
  }

  // Return a copy of the window
  Record window_walker::window()
  {
    return _rec;
  }

  // Skip a number of positions
  void window_walker::skip(size_t n)
  {
    _current_pos += n;
    if(_parent.size() >= (_current_pos + _ws))
      {
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_end = false;
      }
    else if(_include_final && _current_pos < _parent.size())
      {
	_rec = _parent.subseq(_current_pos, _parent.size() - 1);
	_end = false; _include_final = false;
      }
    else
      {
	_end = true;
      }
  }

  // Rewind a number of bases
  void window_walker::rewind(size_t n)
  {
    if(_current_pos < n)
      {
	if(_include_final)
	  {
	    _rec = _parent.subseq(0, _current_pos);
	    _begin = false; _include_final = false;
	  }
	else
	  {
	    _begin = true;
	  }
	_end = _current_pos <= (_parent.size() - 1) ? false : true;
      }
    else
      {
	_current_pos -= n;
	_rec = _parent.subseq(_current_pos, _current_pos + _ws - 1);
	_begin = false;
	_end = _current_pos <= (_parent.size() - 1) ? false : true;
      }
  }
  
};
