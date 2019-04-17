#ifndef VARVAL_BEDPED
#define VARVAL_BEDPE

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/RefGenome.h"

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


struct BEDPEEntry {

  BEDPEEntry(const SeqLib::GenomicRegion& g1, const SeqLib::GenomicRegion& g2) {
    gr1 = g1;
    gr2 = g2;
  }

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

bool BEDPE::ReadBEDPE(const std::string & file, const SeqLib::BamHeader& hdr) {

  gzFile fp = NULL;
  fp = strcmp(file.c_str(), "-")? gzopen(file.c_str(), "r") : gzdopen(fileno(stdin), "r");
  
  if (file.empty() || !fp) {
    std::cerr << "BED file not readable: " << file << std::endl;
    return false;
  }
  
  // http://www.lemoda.net/c/gzfile-read/
  while (1) {
    
    int err;                    
    char buffer[GZBUFFER];
    gzgets(fp, buffer, GZBUFFER);
    int bytes_read = strlen(buffer);

    // get one line
    if (bytes_read < GZBUFFER - 1) {
      if (gzeof (fp)) break;
      else {
	const char * error_string;
	error_string = gzerror (fp, &err);
	if (err) {
	  fprintf (stderr, "Error: %s.\n", error_string);
	  exit (EXIT_FAILURE);
	}
      }
    }
    
    // prepare to loop through each field of BED line
    std::string chr, pos1, pos2, o1;
    std::string Bchr, Bpos1, Bpos2, Bo1;
    std::string line(buffer);
    std::istringstream iss_line(line);
    std::string val;
    if (line.find("#") != std::string::npos) 
      continue;
    
    // read first three BED columns
    iss_line >> chr >> pos1 >> pos2 >> o1 >> Bchr >> Bpos1 >> Bpos2 >> Bo1;

    // construct the GenomicRegion
    SeqLib::GenomicRegion g1(chr, pos1, pos2, hdr);
    SeqLib::GenomicRegion g2(Bchr, Bpos1, Bpos2, hdr);
    g1.strand = o1.at(0);
    g2.strand = Bo1.at(0);

    BEDPEEntry pe(g1, g2);
    pe.id = chr + ":" + pos1 + o1 + "__" + Bchr + ":" + Bpos1 + Bo1;
    
    bedpe.push_back(pe);
  }
  
  return true;
}

std::string BEDPEEntry::convertToContig(SeqLib::RefGenome& ref, const SeqLib::BamHeader& hdr, int pad) const {

  std::string s1, s2;
  if (gr1.strand == '+')
    s1 = ref.QueryRegion(gr1.ChrName(hdr), gr1.pos1 - pad, gr1.pos1);
  else
    s1 = ref.QueryRegion(gr1.ChrName(hdr), gr1.pos1, gr1.pos1 + pad);

  if (gr2.strand == '+')
    s2 = ref.QueryRegion(gr2.ChrName(hdr), gr2.pos1 - pad, gr2.pos1);
  else
    s2 = ref.QueryRegion(gr2.ChrName(hdr), gr2.pos1, gr2.pos1 + pad);

  if (gr1.strand == '-' && gr2.strand == '-')
    return reverseComplement(s2) + s1;
  if (gr1.strand == '+' && gr2.strand == '+')
    return s1 + reverseComplement(s2);
  
  return s1 + s2;

}

std::ostream& operator<<(std::ostream& out, const BEDPE& c) {
  for (auto& i : c.bedpe) 
    out << i << std::endl;
  return out;
}

std::ostream& operator<<(std::ostream& out, const BEDPEEntry& c) {
  
  out << c.gr1 << "----" << c.gr2;
  return out;

}

#endif
