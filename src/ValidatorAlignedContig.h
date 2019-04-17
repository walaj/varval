#ifndef VALIDATOR_ALIGNED_CONTIG_H__
#define VALIDATOR_ALIGNED_CONTIG_H__

#include "SnowTools/AlignedContig.h"

class ValidatorAlignmentFragment : public SnowTools::AlignmentFragment {

 public:
  
 private:

};

class ValidatorAlignedContig : public SnowTools::AlignedContig {

 public: 

 ValidatorAlignedContig(const SnowTools::BamReadVector& bav) : SnowTools::AlignedContig(bav) {}

  int countDiscordants();

  int discordants;

 private:



};



#endif
