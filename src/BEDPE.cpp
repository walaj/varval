#include "BEDPE.h"

BEDPEEntry::BEDPEEntry(const SeqLib::GenomicRegion& g1, const SeqLib::GenomicRegion& g2) {
  gr1 = g1;
  gr2 = g2;
}

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
    iss_line >> chr >> pos1 >> pos2 >> Bchr >> Bpos1 >> Bpos2 >> o1 >> Bo1;
    
    // construct the GenomicRegion
    SeqLib::GenomicRegion g1(chr, pos1, pos2, hdr);
    SeqLib::GenomicRegion g2(Bchr, Bpos1, Bpos2, hdr);
    if (o1.length())
      g1.strand = o1.at(0);
    if (Bo1.length())
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


