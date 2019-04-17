#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <iomanip>


//#include "SnowTools/SVVCF.h"
#include "SeqLib/SeqLibUtils.h"
#include "AlignedContig.h"
#include "BreakPoint.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/SeqLibCommon.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamWalker.h"
#include "SeqLib/BWAWrapper.h"

typedef std::unordered_map<std::string, std::string> BamMap;

std::set<std::string> prefixes;

static const char *VALIDATOR_BAM_USAGE_MESSAGE =
"Usage: validator <input.bam> [OPTIONS] \n\n"
"  Description: Validate variant alleles directly from sequencing reads\n"
"\n"
" General options\n"
"      --help                           Display this help and exit\n"
"  -v, --verbose                        Verbose output\n"
" Required input\n"
"  -f, --vcf-file                       VCF file containing rearrangements\n"
"  -t, --tumor-bam                      Tumor bam file. Can input any number with extra -t calls\n"
"  -n, --normal-bam                     Normal bam file. Can input any number with extra -t calls\n"
"  -G, --reference-genome               Path to indexed reference genome [Broads HG19]\n"
" Optional input\n"
"  -k, --region-file                    Set a region txt file. Format: one region per line, Ex: 1,10000000,11000000\n"
"  -A, --alignment-score-frac           The AS (alignment score) from BWA-MEM must be >= A * read length in order for read to count as aligned to contig. [0.8]\n"
"  -a, --analysis-id                    Set a unique id to prepend output files with. [noid]\n"
"      --express                        \n"
"\n";

std::vector<AlignedContig> alc;
SeqLib::BWAWrapper bwa;
SeqLib::USeqVector usv;

namespace opt {

  static bool verbose = true;
  static std::string vcf;
  static std::string reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

  static BamMap bam;

  static std::string normal_bam;
  static std::string tumor_bam;

  static int sample_number = 0;

  static bool express = false;
  static int pad = 500;

  static std::string analysis_id = "noid";

  static std::string regionFile;

  static double alignment_score_perc = 0.8;

  static double lod = 8;
  static double lod_db = 7;
  static double lod_no_db = 2.5;
  static double lod_germ = 3;

}

enum {
  OPT_HELP
};

static const char* shortopts = "hvG:t:n:f:k:a:eA:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, OPT_HELP },
  { "verbose",                    no_argument, NULL, 'v' },
  { "vcf-file",                   required_argument, NULL, 'f' },
  { "alignment-score-frac",       required_argument, NULL, 'A' },
  { "tumor-bam",                  required_argument, NULL, 't' },
  { "normal-bam",                 required_argument, NULL, 'n' },
  { "region-file",                required_argument, NULL, 'k' },
  { "analysis-id",                required_argument, NULL, 'a' },
  { "express",                    no_argument, NULL, 'e' },
  { "reference-genome",                    required_argument, NULL, 'G' },
  { NULL, 0, NULL, 0 }
};

static struct timespec start;
static ogzstream all_align, all_bps;

// forward declare
void parseVarOptions(int argc, char** argv);
void alignReadsToDerivedContigs(BamWalker& walk, const std::string& prefix);

template <typename T>
  void fopen(const std::string& s, T& o) {
  o.open(s.c_str(), std::ios::out);
}

std::string __bamOptParse(std::unordered_map<std::string, std::string>& obam, std::istringstream& arg, int sample_number, const std::string& prefix) {
  std::stringstream ss;
  std::string bam;
    arg >> bam;
    ss.str(std::string());
    ss << prefix << std::setw(3) << std::setfill('0') << sample_number;
    obam[ss.str()] = bam;
    prefixes.insert(ss.str());
    return bam;
}

int main(int argc, char** argv) {

#ifndef __APPLE__
  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
#endif

  // parse the command line
  parseVarOptions(argc, argv);

  fopen(opt::analysis_id + ".alignments.txt.gz", all_align);
  fopen(opt::analysis_id + ".bps.txt.gz", all_bps);

  all_bps << BreakPoint::validator_header();


  // start the timer
  clock_gettime(CLOCK_MONOTONIC, &start);
  
  // read a tumor BAM
  BamWalker tum, norm;
  if (!opt::normal_bam.empty()) {
    std::cerr << "...opening normal BAM " << opt::normal_bam << std::endl;
    norm =  BamWalker(opt::normal_bam);
  }
  if (!opt::tumor_bam.empty()) {
    std::cerr << "...opening tumor BAM " << opt::tumor_bam << std::endl;
    tum = SeqLib::BamReader(opt::tumor_bam);
  }
  
  // add the bps headers
  for (auto& b : opt::bam) 
    all_bps << "\t" << b.first << "_" << boost::filesystem::path(opt::tumor_bam).filename().string();

  all_bps << std::endl;

  // parse the region file
  GRC regions;
  bool rgfile = true; //SnowTools::read_access_test(opt::regionFile);
    
  int pad = opt::express ? 1000 : 0;

  if (rgfile)
    regions.regionFileToGRV(opt::regionFile, pad);
  // parse as a samtools string eg 1:1,000,000-2,000,000
  else if (opt::regionFile.find(":") != std::string::npos && opt::regionFile.find("-") != std::string::npos)
    regions.add(SnowTools::GenomicRegion(opt::regionFile, tum.header()));
  
  // it's a single chromosome
  else if (!opt::regionFile.empty()) {
    SnowTools::GenomicRegion gr(opt::regionFile, "1", "1", tum.header());
    if (gr.chr == -1 || gr.chr >= tum.header()->n_targets) {
      std::cerr << "ERROR: Trying to match chromosome " << opt::regionFile << " to one in header, but not match found" << std::endl;
      exit(EXIT_FAILURE);
    } else {
      gr.pos2 = tum.header()->target_len[gr.chr];
      regions.add(gr);
    }
  }

  if (regions.size()) {
    if (!opt::tumor_bam.empty())
      tum.setBamWalkerRegions(regions.asGenomicRegionVector());
    if (!opt::normal_bam.empty())
      norm.setBamWalkerRegions(regions.asGenomicRegionVector());
  }

  SVVCF svvcf(opt::vcf, tum.header());

  // read an index
  RefGenome ref(opt::reference);

  if (opt::verbose)
    std::cerr << "...read VCF file " << opt::vcf << " with " << svvcf.size() << " entries" << std::endl;

   // loop through regions and make split references}
  for (size_t i = 0; i < svvcf.size(); ++i) {
    for (size_t rc_ins = 0; rc_ins < 2; ++rc_ins) {
      
      BamReadVector brv;
      std::string name = svvcf.getVariantID(i) + "_INSRC_" + std::to_string(rc_ins); 
      std::string seq = svvcf.getSurroundingContigFromReference(i, ref, tum.header(), brv, opt::pad, rc_ins == 1);
      if (seq.empty())
	continue;
      usv.push_back({name, seq}); 
      
      SnowTools2::AlignedContig ac(brv, prefixes);
      ac.forceBreakPoint(svvcf.getSide1(i), svvcf.getSide2(i));
      alc.push_back(ac);
    }
  }
  
  // construct the new reference
  bwa.constructIndex(usv);
  
  std::vector<int> count_vec(svvcf.size(), 0);

  for (auto& b : opt::bam) 
    if (b.first.at(0) == 't')
      alignReadsToDerivedContigs(tum, b.first);
    else 
      alignReadsToDerivedContigs(norm, b.first);
  
  for (auto& a : alc) {
    a.decoyReads(tum.header(), ref);
    a.countDiscordants();
    a.assignSupportCoverage();
    a.splitCoverage();	
  }

  size_t cc = 0;
  for (auto& a : alc)
    cc += a.m_bamreads.size();

  std::cerr << "...counted total of "  << cc << " reads across " << alc.size() << " SV" << std::endl;
  std::cerr << SnowTools::displayRuntime(start) << std::endl;

  std::vector<SnowTools2::BreakPoint> bp_glob, bp_glob_new;
  for (auto& a : alc) {
    std::vector<SnowTools2::BreakPoint> allbreaks = a.getAllBreakPoints();
    bp_glob.insert(bp_glob.end(), allbreaks.begin(), allbreaks.end());
  }

  // dedube
  /*
  std::unordered_map<std::string, int> alt_max;
  for(auto& b : bp_glob)
    alt_max[b.cname] = 0;

  for (auto& b : bp_glob) {
    for (auto& a : b.allele)
      alt_max[b.cname] = std::max(alt_max[b.cname], a.second.alt);
  }
    
  for (auto& b : bp_glob) {
    size_t alt_max_this = 0;
    for (auto& a : b.allele)    
      alt_max_this = std::max(alt_max_this, a.second.alt);
    if (alt_max[b.cname] == alt_max_this)
      bp_glob_new.push_back(bp.glob);
  }
  */
  std::cerr << " bp.glob.size() " << bp_glob.size() << std::endl;

  // set homologies
  for (auto& b : bp_glob)
    b.__set_homologies_insertions();

  // repeat sequence filter
  for (auto& i : bp_glob)
    i.repeatFilter();

  for (auto& i : bp_glob)
    i.scoreBreakpoint(opt::lod, opt::lod_db, opt::lod_no_db, opt::lod_germ);

  bool PRINT_READS = true;

  // de duplicate the breakpoints
  std::sort(bp_glob.begin(), bp_glob.end());
  bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );
  
  for (auto& i : bp_glob) {
    i.scrubName();
    all_bps << i.toValidatorFileString(!PRINT_READS) << std::endl;
  }

  for (auto& a : alc) {
    a.scrubName();
    all_align << a << std::endl;
  }
  
  // write the reads to a BAM
  SnowTools::BamWalker nbam, tbam;

  if (!opt::tumor_bam.empty()) {
    bam_hdr_t * r2c_hdr = bam_hdr_dup(tum.header());                                                                                                                                                                                                                                                                                                                                      
    tbam.SetWriteHeader(r2c_hdr);                                                                                                                                                                                                                                                                                                                                                             
    tbam.OpenWriteBam(opt::analysis_id + ".tumor.r2contig.bam");
  }

  if (!opt::normal_bam.empty()) {
    bam_hdr_t * r2c_hdrn = bam_hdr_dup(norm.header());                                                                                                                                                                                                                                                                                                                                      
    nbam.SetWriteHeader(r2c_hdrn);                                                                                                                                                                                                                                                                                                                                                             
    nbam.OpenWriteBam(opt::analysis_id + ".normal.r2contig.bam");
  }

  std::unordered_set<std::string> seen;
  for (auto& b : bp_glob) {
    for (auto& r : b.reads) {
      std::string sr = r.GetZTag("SR");
      if (!seen.count(sr)) { 
	seen.insert(sr);
	r.AddZTag("CC", b.cname + "C");
	if (sr.at(0) == 't')
	  tbam.writeAlignment(r);
	else
	  nbam.writeAlignment(r);

      }
    }
  }
    

}

void parseVarOptions(int argc, char** argv) {
  
  bool die = false;

  if (argc < 2) {
    std::cerr << "\n" << VALIDATOR_BAM_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {

    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case OPT_HELP: die = true; break;
    case 'v': opt::verbose = true; break;
    case 'f': arg >> opt::vcf; break;
    case 'e': opt::express = true; break;
    case 'G': arg >> opt::reference; break;
    case 'A': arg >> opt::alignment_score_perc; break;
    case 'k': arg >> opt::regionFile; break;
    case 'a': arg >> opt::analysis_id; break;
    case 't': 
      tmp = __bamOptParse(opt::bam, arg, opt::sample_number++, "t");
      if (opt::tumor_bam.empty())
	opt::tumor_bam = tmp;
      if (!read_access_test(tmp)) {
	std::cerr << "Tumor BAM file " << tmp << " does not exist / not readable " << std::endl;
	die = true;
      }
      break;
    case 'n': 
      tmp = __bamOptParse(opt::bam, arg, opt::sample_number++, "n");
      if (opt::normal_bam.empty())
	opt::normal_bam = tmp;
      if (!read_access_test(tmp)) {
	std::cerr << "Normal BAM file " << tmp << " does not exist / not readable " << std::endl;
	die = true;
      }
      break;
    }
  }

  if (!opt::vcf.empty() && !read_access_test(opt::vcf)) {
    std::cerr << "VCF file: " << opt::vcf << " is not available" << std::endl;
    die = true;
  }

  if (!read_access_test(opt::reference)) {
    std::cerr << "Reference genome file: " << opt::reference << " does not exist " << std::endl;
    die = true;
  }

  
  if (opt::express)
    opt::regionFile = opt::vcf;
  
  // something went wrong, kill
  if (die) {
    std::cerr << "\n" << VALIDATOR_BAM_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

}

void alignReadsToDerivedContigs(BamWalker& walk, const std::string& prefix) {

  // 
  BamRead i; 
  bool rule;
  size_t count = 0;
  if (opt::verbose)
    std::cerr << "...looping through the " << prefix << " BAM" << std::endl;

  while (walk.GetNextRead(i, rule)) {
    
    // make sure we have only one alignment per contig
    std::set<std::string> cc;

    ++count;
    if (count % 100000 == 0) {
      size_t cc = 0;
      for (auto& a : alc)
	cc += a.m_bamreads.size();
      std::cerr << "... at read " << AddCommas(count) << " at pos " << i.Brief(walk.header()) << " aligned " << AddCommas(cc) << " read ";
      std::cerr << SnowTools::displayRuntime(start) << std::endl;
    }

    BamReadVector brv;
    std::string seqr = i.QualitySequence();
    bwa.alignSingleSequence(seqr, i.Qname(), brv, true, 0.90, 50); // hardclip = true

    if (brv.size() == 0) 
      continue;

    // add the ID tag
    std::string srn = prefix + "_" + std::to_string(i.AlignmentFlag()) + "_" + i.Qname();
    assert(srn.length());
    i.AddZTag("SR", srn);

    BamReadVector bpass;

    for (auto& r : brv) { // r is read aligned to contigs
      
      bool length_pass = (r.PositionEnd() - r.Position()) >= ((double)seqr.length() * 0.75);
      
      if (length_pass && !cc.count(usv[r.ChrID()].name) && r.NumClip() < 10 && (double)r.GetIntTag("AS") > opt::alignment_score_perc * (double)seqr.length()) {
	bpass.push_back(r);
	cc.insert(usv[r.ChrID()].name);
      }
    }

    // annotate the original read
    for (auto& r : bpass) {

      if (r.ReverseFlag())
	i.SmartAddTag("RC","1");
      else 
	i.SmartAddTag("RC","0");

      i.SmartAddTag("HC", r.Sequence()); // add the hardclip
      i.SmartAddTag("AA", std::to_string(r.GetIntTag("AS"))); // add the alignment score
      i.SmartAddTag("SL", std::to_string(r.Position()));
      i.SmartAddTag("SE", std::to_string(r.PositionEnd()));
      i.SmartAddTag("TS", std::to_string(r.AlignmentPosition()));
      i.SmartAddTag("TE", std::to_string(r.AlignmentEndPosition()));
      i.SmartAddTag("SC", r.CigarString());
      i.SmartAddTag("CN", usv[r.ChrID()].name/*getContigName()*/);
      
      for (auto& a : alc) { // assign to a contig

	if (a.getContigName() != usv[r.ChrID()].name) 
	  continue;

	a.m_bamreads.push_back(i);
	
	// add the coverage to the aligned contig
	int cc = r.Position();
	std::string srr = i.GetZTag("SR");
	while (cc <= r.PositionEnd() && cc < (int)a.tum_cov.size()) 
	  ++a.cov[srr.substr(0,4)][cc++];
	if (srr.at(0) == 't') 
	  while (cc <= r.PositionEnd() && cc < (int)a.tum_cov.size()) 
	    ++a.tum_cov[cc++];
	else 
	  while (cc <= r.PositionEnd() && cc < (int)a.norm_cov.size())
	    ++a.norm_cov[cc++];
      } // alc lllop
      
    } // bpas loop
    
  }


}

