#ifndef VARVAL_BEDPED_H
#define VARVAL_BEDPED_H

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/RefGenome.h"

/*
char complement(char n)
{   
  switch(n)
    {   
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    }   
  assert(false);
  return ' ';
}   

// Complement a sequence
std::string complement(const std::string& seq)
{
  std::string out(seq.length(), 'A');
  size_t l = seq.length();
  for(size_t i = 0; i < l; ++i)
    out[i] = complement(seq[i]);
  return out;
}

// Reverse complement a sequence
std::string reverseComplement(const std::string& seq)
{
  std::string out(seq.length(), 'A');
  size_t last_pos = seq.length() - 1;
  for(int i = last_pos; i >= 0; --i)
    {
      out[last_pos - i] = complement(seq[i]);
    }
  return out;
}
*/

struct BEDPEEntry {

  BEDPEEntry(const SeqLib::GenomicRegion& g1, const SeqLib::GenomicRegion& g2);

  SeqLib::GenomicRegion gr1;
  SeqLib::GenomicRegion gr2;

  std::string id;

  std::string convertToContig(SeqLib::RefGenome& ref, const SeqLib::BamHeader& hdr, int pad) const;

  friend std::ostream& operator<<(std::ostream& out, const BEDPEEntry& c);

};

struct BEDPE {

  bool ReadBEDPE(const std::string & file, const SeqLib::BamHeader& hdr);

  std::vector<BEDPEEntry> bedpe;

  friend std::ostream& operator<<(std::ostream& out, const BEDPE& c);
};

#endif
