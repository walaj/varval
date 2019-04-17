#include "ValidatorAlignedContig.h"

int ValidatorAlignedContig::countDiscordants() {

  std::pair<int,int> break_pair_1 = m_frag_v[0].getContigBreaks();
  std::pair<int,int> break_pair_2 = m_frag_v[1].getContigBreaks();

  int cpos1 = break_pair_1.second;
  int cpos2 = break_pair_2.first;

  std::unordered_map<std::string, std::pair<int,int>> all;

  for (auto& r : m_bamreads) {

      // get the more complex tags (since there can be multiple annotations per tag)
    std::vector<int> posvec = r.GetSmartIntTag("SL"); // start positions ON CONTIG
    std::vector<std::string> cnvec = r.GetSmartStringTag("CN");
    std::vector<int> alnvec = r.GetSmartIntTag("TS"); // start positions ON READ
    std::vector<std::string> hcvec = r.GetSmartStringTag("HC"); // hard-clipped seq on CONTIG      
    std::vector<int> rcvec = r.GetSmartIntTag("RC"); // read reverse complmented relative to contig
    if (posvec.size() != cnvec.size() || alnvec.size() != posvec.size() || hcvec.size() != posvec.size())
      continue;
    int pos = 0, aln = 0, rc;
    std::string hc;
    size_t kk = 0;
    for (; kk < cnvec.size(); kk++) 
      if (cnvec[kk] == getContigName()) {
	pos = posvec[kk];
	aln = alnvec[kk];
	rc = rcvec[kk];
	if (hcvec.size())
	  hc = hcvec[kk];
	break;
      }
    
    if (pos + 30  < cpos1) { // end part fits?
      std::string qq = r.Qname();
      int flag = r.FirstFlag() ? 1 : 2;
      assert(flag == 1 || flag == 2);
      all[qq].first = flag;
    } else if (pos > cpos2 - 5) {
      std::string qq = r.Qname();
      int flag = r.FirstFlag() ? 1 : 2;      
      assert(flag == 1 || flag == 2);
      all[qq].second = flag;
    }
  }

  discordants = 0;
  for(auto& i : all)
    if (i.second.first != i.second.second && i.second.first + i.second.second == 3)
      ++discordants;
  return discordants;
}

