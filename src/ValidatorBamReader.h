#ifndef VALIDATOR_BAM_WALKER_H__
#define VALIDATOR_BAM_WALKER_H__

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "svabaRead.h"

class validatorBamReader: public SeqLib::BamReader {

 public:
  
  validatorBamReader() {}

  bool GetNextFilteredRecord(svabaRead& s);

  void readBam();

  void QualityTrimRead(svabaRead& r) const;

  std::string prefix;

};
  

#endif
