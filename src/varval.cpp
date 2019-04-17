#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <set>
#include "gzstream.h"

#include "BEDPE.h"
#include "SeqLib/BamReader.h"

static const char *VARVAL_BAM_USAGE_MESSAGE =
"Usage: varval <input.bam> [OPTIONS] \n\n"
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

typedef std::unordered_map<std::string, std::string> BamMap;

std::set<std::string> prefixes;

static struct timespec start;
static ogzstream all_align, all_bps;

namespace opt {

  static bool verbose = true;
  static std::string bedpe;
  static std::string reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

  static SeqLib::BamHeader hdr;

  static BamMap bam;

  static int sample_number = 0;

  static bool express = false;
  static int pad = 500;

  static SeqLib::RefGenome ref;

  static std::string analysis_id = "noid";

  static std::string regionFile;

  static std::string normal_bam;
  static std::string tumor_bam;

}

enum {
  OPT_HELP
};

static const char* shortopts = "hvG:t:n:f:k:a:eA:";
static const struct option longopts[] = {
  { "help",                       no_argument, NULL, OPT_HELP },
  { "verbose",                    no_argument, NULL, 'v' },
  { "bedpe-file",                   required_argument, NULL, 'f' },
  { "alignment-score-frac",       required_argument, NULL, 'A' },
  { "tumor-bam",                  required_argument, NULL, 't' },
  { "normal-bam",                 required_argument, NULL, 'n' },
  { "region-file",                required_argument, NULL, 'k' },
  { "analysis-id",                required_argument, NULL, 'a' },
  { "express",                    no_argument, NULL, 'e' },
  { "reference-genome",                    required_argument, NULL, 'G' },
  { NULL, 0, NULL, 0 }
};


// forward declare
void parseVarOptions(int argc, char** argv);

template <typename T>
  void fopen(const std::string& s, T& o) {
  o.open(s.c_str(), std::ios::out);
}

std::string __bamOptParse(BamMap& obam, std::istringstream& arg, 
        int sample_number, const std::string& prefix) {

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

  // open reference genome
  opt::ref.LoadIndex(opt::reference);

  // open the BED file
  std::cerr << "parsing input " << opt::bedpe << std::endl;
  BEDPE pe;
  pe.ReadBEDPE(opt::bedpe, opt::hdr);

  for (auto& i : pe.bedpe) {
    std::cout << ">" << i.id << std::endl;
    std::cout << i.convertToContig(opt::ref, opt::hdr, opt::pad) << std::endl;
  }

  return 0;

}

void parseVarOptions(int argc, char** argv) {
  
  bool die = false;

  SeqLib::BamReader br;

  if (argc < 2) {
    std::cerr << "\n" << VARVAL_BAM_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {

    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case OPT_HELP: die = true; break;
    case 'v': opt::verbose = true; break;
    case 'f': arg >> opt::bedpe; break;
    case 'e': opt::express = true; break;
    case 'G': arg >> opt::reference; break;
    case 'k': arg >> opt::regionFile; break;
    case 'a': arg >> opt::analysis_id; break;
    case 't': 
      tmp = __bamOptParse(opt::bam, arg, opt::sample_number++, "t");
      if (opt::tumor_bam.empty())
	opt::tumor_bam = tmp;
      br.Open(tmp);
      opt::hdr = br.Header();
      //if (!read_access_test(tmp)) {
      //	std::cerr << "Tumor BAM file " << tmp << " does not exist / not readable " << std::endl;
      //	die = true;
      //}
      break;
    case 'n': 
      tmp = __bamOptParse(opt::bam, arg, opt::sample_number++, "n");
      if (opt::normal_bam.empty())
	opt::normal_bam = tmp;
      //if (!read_access_test(tmp)) {
      //std::cerr << "Normal BAM file " << tmp << " does not exist / not readable " << std::endl;
	//die = true;
      //}
      break;
    }
  }

  //if (!opt::bedpe.empty()) { // && !read_access_test(opt::vcf)) {
  // std::cerr << "BEDPE file: " << opt::bedpe << " is not available" << std::endl;
  // die = true;
  //}

  //if (!read_access_test(opt::reference)) {
  //  std::cerr << "Reference genome file: " << opt::reference << " does not exist " << std::endl;
  //  die = true;
  // }

  
  if (opt::express)
    opt::regionFile = opt::bedpe;
  
  // something went wrong, kill
  if (die) {
    std::cerr << "\n" << VARVAL_BAM_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

}
