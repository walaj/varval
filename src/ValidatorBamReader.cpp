#include "ValidatorBamReader.h"

bool validatorBamReader::GetNextFilteredRecord(svabaRead& s) {

  SeqLib::BamRecord r;
  if (!GetNextRecord(r))
    return false;
  
  s = svabaRead(r, prefix);

  // quality score trim read
  QualityTrimRead(s);

  return true;
}

void validatorBamReader::QualityTrimRead(svabaRead& r) const {

  int32_t startpoint = 0, endpoint = 0;
  r.QualityTrimmedSequence(3, startpoint, endpoint);
  int32_t new_len = endpoint - startpoint;

  // read is not good
  if (endpoint != -1 && new_len < r.Length() && new_len > 0 && new_len - startpoint >= 0 && startpoint + new_len <= r.Length()) { 
    try { 
      r.SetSeq(r.Sequence().substr(startpoint, new_len));
    } catch (...) {
      std::cerr << "Subsequence failure with sequence of length "  
		<< r.Sequence().length() << " and startpoint "
		<< startpoint << " endpoint " << endpoint 
		<< " newlen " << new_len << std::endl;
    }

  // read is fine
  } else {
    r.SetSeq(r.Sequence()); // copies the sequence
  }
  

}
