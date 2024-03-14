//
// Created by hlf on 6/13/22.
//

#ifndef MUTECT2CPP_MASTER_VCFWRITER_H
#define MUTECT2CPP_MASTER_VCFWRITER_H

#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "VariantContext.h"
#include "samtools/SAMSequenceDictionary.h"
#include "VCFConstants.h"

#include "annotator/BaseQuality.h"
#include "annotator/FragmentLength.h"
#include "annotator/Coverage.h"
#include "annotator/MappingQuality.h"
#include "annotator/ReadPosition.h"

#include "annotator/DepthPerSampleHC.h"
#include "annotator/DepthPerAlleleBySample.h"
#include "annotator/OrientationBiasReadCounts.h"
#include "annotator/StrandBiasBySample.h"


class VCFWriter {
    phmap::flat_hash_map<std::string, std::string> mInfoMetaData = {{VCFConstants::MEDIAN_BASE_QUALITY_KEY,         "INFO=<ID=MBQ,Number=R,Type=Integer,Description=\"median base quality\">"},
	                                                              {VCFConstants::DEPTH_KEY,                       "INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">"},
	                                                              {VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY,      "INFO=<ID=MFRL,Number=R,Type=Integer,Description=\"median fragment length\">"},
	                                                              {VCFConstants::MEDIAN_MAPPING_QUALITY_KEY,      "INFO=<ID=MMQ,Number=R,Type=Integer,Description=\"median mapping quality\">"},
	                                                              {VCFConstants::MEDIAN_READ_POSITON_KEY,         "INFO=<ID=MPOS,Number=A,Type=Integer,Description=\"median distance from end of read\">"},
	                                                              {VCFConstants::STR_PRESENT_KEY,                 "INFO=<ID=STR,Number=0,Type=Flag,Description=\"Variant is a short tandem repeat\">"},
	                                                              {VCFConstants::REPEAT_UNIT_KEY,                 "INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">"},
	                                                              {VCFConstants::NORMAL_LOG_10_ODDS_KEY,          "INFO=<ID=NLOD,Number=A,Type=Float,Description=\"Normal log 10 likelihood ratio of diploid het or hom alt genotypes\">"},
	                                                              {VCFConstants::REPEATS_PER_ALLELE_KEY,          "INFO=<ID=RPA,Number=.,Type=Integer,Description=\"Number of times tandem repeat unit is repeated, for each allele (including reference)\">"},
	                                                              {VCFConstants::TUMOR_LOG_10_ODDS_KEY,           "INFO=<ID=TLOD,Number=A,Type=Float,Description=\"Log 10 likelihood ratio score of variant existing versus not existing\">"},
	                                                              {VCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY, "INFO=<ID=NALOD,Number=A,Type=Float,Description=\"Negative log 10 odds of artifact in normal with same allele fraction as tumor\">"},
	                                                              {VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY,    "INFO=<ID=ECNT,Number=1,Type=Integer,Description=\"Number of events in this haplotype\">"},
	                                                              {VCFConstants::IN_PON_KEY,                      "INFO=<ID=PON,Number=0,Type=Flag,Description=\"site found in panel of normals\">"},
	                                                              {VCFConstants::POPULATION_AF_KEY,               "INFO=<ID=POPAF,Number=A,Type=Float,Description=\"negative log 10 population allele frequencies of alt alleles\">"},
	                                                              {VCFConstants::GERMLINE_QUAL_KEY,               "INFO=<ID=GERMQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not germline variants\">"},
	                                                              {VCFConstants::CONTAMINATION_QUAL_KEY,          "INFO=<ID=CONTQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to contamination\">"},
	                                                              {VCFConstants::SEQUENCING_QUAL_KEY,             "INFO=<ID=SEQQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not sequencing errors\">"},
	                                                              {VCFConstants::POLYMERASE_SLIPPAGE_QUAL_KEY,    "INFO=<ID=STRQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors\">"},
	                                                              {VCFConstants::READ_ORIENTATION_QUAL_KEY,       "INFO=<ID=ROQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to read orientation artifact\">"},
	                                                              {VCFConstants::STRAND_QUAL_KEY,                 "INFO=<ID=STRANDQ,Number=1,Type=Integer,Description=\"Phred-scaled quality of strand bias artifact\">"},
	                                                              {VCFConstants::ORIGINAL_CONTIG_MISMATCH_KEY,    "INFO=<ID=OCM,Number=1,Type=Integer,Description=\"Number of alt reads whose original alignment doesn't match the current contig.\">"},
	                                                              {VCFConstants::N_COUNT_KEY,                     "INFO=<ID=NCount,Number=1,Type=Integer,Description=\"Count of N bases in the pileup\">"},
	                                                              {VCFConstants::UNIQUE_ALT_READ_SET_COUNT_KEY,   "INFO=<ID=UNIQ_ALT_READ_COUNT,Number=1,Type=Integer,Description=\"Number of ALT reads with unique start and mate end positions at a variant site\">"}};

    phmap::flat_hash_map<std::string, std::string> mFormatMetaData = {{VCFConstants::GENOTYPE_ALLELE_DEPTHS,          "FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"},
	                                                                {VCFConstants::DEPTH_KEY,                       "FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">"},
	                                                                {VCFConstants::F1R2_KEY,                        "FORMAT=<ID=F1R2,Number=R,Type=Integer,Description=\"Count of reads in F1R2 pair orientation supporting each allele\">"},
	                                                                {VCFConstants::F2R1_KEY,                        "FORMAT=<ID=F2R1,Number=R,Type=Integer,Description=\"Count of reads in F2R1 pair orientation supporting each allele\">"},
	                                                                {VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY,       "FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">"},
	                                                                {VCFConstants::GENOTYPE_KEY,                    "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"},
	                                                                {VCFConstants::GENOTYPE_QUALITY_KEY,            "FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"},
	                                                                {VCFConstants::GENOTYPE_PL_KEY,                 "FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"},
	                                                                {VCFConstants::ALLELE_FRACTION_KEY,             "FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">"},
	                                                                {VCFConstants::HAPLOTYPE_CALLER_PHASING_ID_KEY, "FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">"},
	                                                                {VCFConstants::HAPLOTYPE_CALLER_PHASING_GT_KEY, "FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles\">"},
	                                                                {VCFConstants::PHASE_SET_KEY,                   "FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing set (typically the position of the first variant in the set)\">"}};

    phmap::flat_hash_set<std::string> mFilterMetaData = {"##FILTER=<ID=base_qual,Description=\"alt median base quality\">",
    "##FILTER=<ID=clustered_events,Description=\"Clustered events observed in the tumor\">",
    "##FILTER=<ID=contamination,Description=\"contamination\">",
    "##FILTER=<ID=duplicate,Description=\"evidence for alt allele is overrepresented by apparent duplicates\">",
    "##FILTER=<ID=fragment,Description=\"abs(ref - alt) median fragment length\">",
    "##FILTER=<ID=germline,Description=\"Evidence indicates this site is germline, not somatic\">",
    "##FILTER=<ID=haplotype,Description=\"Variant near filtered variant on same haplotype.\">",
    "##FILTER=<ID=low_allele_frac,Description=\"Allele fraction is below specified threshold\">",
    "##FILTER=<ID=map_qual,Description=\"ref - alt median mapping quality\">",
    "##FILTER=<ID=multiallelic,Description=\"Site filtered because too many alt alleles pass tumor LOD\">",
    "##FILTER=<ID=n_ratio,Description=\"Ratio of N to alt exceeds specified ratio\">",
    "##FILTER=<ID=normal_artifact,Description=\"artifact_in_normal\">",
    "##FILTER=<ID=numt_chimera,Description=\"NuMT variant with too many ALT reads originally from autosome\">",
    "##FILTER=<ID=numt_novel,Description=\"Alt depth is below expected coverage of NuMT in autosome\">",
    "##FILTER=<ID=orientation,Description=\"orientation bias detected by the orientation bias mixture model\">",
    "##FILTER=<ID=panel_of_normals,Description=\"Blacklisted site in panel of normals\">",
    "##FILTER=<ID=position,Description=\"median distance of alt variants from end of reads\">",
    "##FILTER=<ID=slippage,Description=\"Site filtered due to contraction of short tandem repeat region\">",
    "##FILTER=<ID=strand_bias,Description=\"Evidence for alt allele comes from one read direction only\">",
    "##FILTER=<ID=strict_strand,Description=\"Evidence for alt allele is not represented in both directions\">",
    "##FILTER=<ID=weak_evidence,Description=\"Mutation does not meet likelihood threshold\">"};

    std::vector<std::string> index2filter = {"PASS", "low_allele_frac", "panel_of_normals", "strand_bias", "germline", "n_ratio", "strict_strand",
                                             "numt_chimera", "slippage", "weak_evidence", "position", "contamination", "multiallelic", "fragment",
                                             "base_qual", "haplotype", "clustered_events", "numt_novel", "duplicate", "map_qual", "orientation", "normal_artifact"};

	std::vector<std::string> contigMetaData;
    phmap::flat_hash_map<std::string, std::string> mOtherMetaData;
	std::vector<std::string> sampleNamesInOrder;
	std::vector<std::string> metaDataInSortedOrder;
	std::string outPath;
	SAMSequenceDictionary refDict;
	htsFile *outFile;
	bcf_hdr_t *hdr;

public:

	VCFWriter(const std::string &outPath, SAMSequenceDictionary sequenceDictionary);

	static std::string getCommandLine(int argc, char *argv[]);

	void
	writeHeader(const std::string &cmdLine, const std::vector<SAMReadGroupRecord> &readGroup, const std::string& normalSample);

	void add(std::shared_ptr<VariantContext>& vc);

	void close();

private:

	void appendHeader_c(char str[]);

	void appendHeader(const std::string &str);
};

#endif //MUTECT2CPP_MASTER_VCFWRITER_H
