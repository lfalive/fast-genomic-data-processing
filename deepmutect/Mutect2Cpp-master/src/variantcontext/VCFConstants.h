//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_VCFCONSTANTS_H
#define MUTECT2CPP_MASTER_VCFCONSTANTS_H

#include <string>

class VCFConstants {
public:
	// reserved INFO/FORMAT field keys
	inline const static std::string ANCESTRAL_ALLELE_KEY = "AA";    //---in c++17, you can use inline const static like Java
	inline const static std::string ALLELE_COUNT_KEY = "AC";
	inline const static std::string ALLELE_FREQUENCY_KEY = "AF";
	inline const static std::string ALLELE_NUMBER_KEY = "AN";
	inline const static std::string RMS_BASE_QUALITY_KEY = "BQ";
	inline const static std::string CIGAR_KEY = "CIGAR";
	inline const static std::string DBSNP_KEY = "DB";
	inline const static std::string DEPTH_KEY = "DP";
	inline const static std::string END_KEY = "END";

	inline const static std::string GENOTYPE_FILTER_KEY = "FT";
	inline const static std::string GENOTYPE_KEY = "GT";
	inline const static std::string GENOTYPE_POSTERIORS_KEY = "GP";
	inline const static std::string GENOTYPE_QUALITY_KEY = "GQ";
	inline const static std::string GENOTYPE_ALLELE_DEPTHS = "AD"; //AD isn't reserved, but is specifically handled by VariantContext
	inline const static std::string GENOTYPE_PL_KEY = "PL";   // phred-scaled genotype likelihoods
	inline const static std::string EXPECTED_ALLELE_COUNT_KEY = "EC";
	inline const static std::string PHASE_SET_KEY = "PS";
	inline const static std::string STR_PRESENT_KEY = "STR";
	inline const static std::string REFSAMPLE_DEPTH_KEY = "REFDEPTH";
	inline const static std::string REPEAT_UNIT_KEY = "RU";
	inline const static std::string REPEATS_PER_ALLELE_KEY = "RPA";

	// missing/default values
	inline const static std::string UNFILTERED = ".";
	inline const static std::string PASSES_FILTERS_v3 = "0";
	inline const static std::string PASSES_FILTERS_v4 = "PASS";
	inline const static std::string EMPTY_ID_FIELD = ".";
	inline const static std::string EMPTY_INFO_FIELD = ".";
	inline const static std::string EMPTY_ALTERNATE_ALLELE_FIELD = ".";
	inline const static std::string MISSING_VALUE_v4 = ".";
	inline const static std::string MISSING_QUALITY_v3 = "-1";
	inline const static double MISSING_QUALITY_v3_DOUBLE = stod(MISSING_QUALITY_v3);

	// Mutect2-specific INFO keys
	inline const static std::string TUMOR_LOG_10_ODDS_KEY = "TLOD";
	inline const static std::string POPULATION_AF_KEY = "POPAF";
	inline const static std::string NORMAL_LOG_10_ODDS_KEY = "NLOD";
	inline const static std::string NORMAL_ARTIFACT_LOG_10_ODDS_KEY = "NALOD";
	inline const static std::string MEDIAN_BASE_QUALITY_KEY = "MBQ";
	inline const static std::string MEDIAN_MAPPING_QUALITY_KEY = "MMQ";
	inline const static std::string MEDIAN_FRAGMENT_LENGTH_KEY = "MFRL";
	inline const static std::string MEDIAN_READ_POSITON_KEY = "MPOS";
	inline const static std::string IN_PON_KEY = "PON";
	inline const static std::string GERMLINE_QUAL_KEY = "GERMQ";
	inline const static std::string CONTAMINATION_QUAL_KEY = "CONTQ";
	inline const static std::string SEQUENCING_QUAL_KEY = "SEQQ";
	inline const static std::string POLYMERASE_SLIPPAGE_QUAL_KEY = "STRQ";
	inline const static std::string STRAND_QUAL_KEY = "STRANDQ";
	inline const static std::string READ_ORIENTATION_QUAL_KEY = "ROQ";
	inline const static std::string ORIGINAL_CONTIG_MISMATCH_KEY = "OCM";
	inline const static std::string N_COUNT_KEY = "NCount";
	inline const static std::string UNIQUE_ALT_READ_SET_COUNT_KEY = "UNIQ_ALT_READ_COUNT";

	// FORMAT keys
	inline const static std::string HAPLOTYPE_CALLER_PHASING_ID_KEY = "PID";
	inline const static std::string HAPLOTYPE_CALLER_PHASING_GT_KEY = "PGT";
	inline const static std::string STRAND_BIAS_BY_SAMPLE_KEY = "SB";

	// M2-specific FORMAT keys
	inline const static std::string ALLELE_FRACTION_KEY = "AF";

	// INFO keys
	inline const static std::string EVENT_COUNT_IN_HAPLOTYPE_KEY = "ECNT"; //M2
	inline const static std::string F1R2_KEY = "F1R2";
	inline const static std::string F2R1_KEY = "F2R1";

    inline const static std::string BAD_HAPLOTYPE_FILTER_NAME = "haplotype";
    inline const static std::string TUMOR_EVIDENCE_FILTER_NAME = "weak_evidence";
    inline const static std::string STRAND_ARTIFACT_FILTER_NAME = "strand_bias";
    inline const static std::string MEDIAN_BASE_QUALITY_FILTER_NAME = "base_qual";
    inline const static std::string MEDIAN_MAPPING_QUALITY_FILTER_NAME = "map_qual";
    inline const static std::string DUPLICATED_EVIDENCE_FILTER_NAME = "duplicate";
    inline const static std::string ARTIFACT_IN_NORMAL_FILTER_NAME = "normal_artifact";
    inline const static std::string GERMLINE_RISK_FILTER_NAME = "germline";
    inline const static std::string CLUSTERED_EVENTS_FILTER_NAME = "clustered_events";
    inline const static std::string LOW_QUAL_FILTER_NAME = "LowQual";
    inline const static std::string ALIGNMENT_ARTIFACT_FILTER_NAME = "alignment";
    inline const static std::string PON_FILTER_NAME = "panel_of_normals";
    inline const static std::string N_RATIO_FILTER_NAME = "n_ratio";
    inline const static std::string READ_POSITION_FILTER_NAME = "position";
    inline const static std::string ALLELE_FRACTION_FILTER_NAME = "low_allele_frac";
    inline const static std::string MULTIALLELIC_FILTER_NAME = "multiallelic";
    inline const static std::string MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME = "fragment";
    inline const static std::string POLYMERASE_SLIPPAGE = "slippage";
};


#endif //MUTECT2CPP_MASTER_VCFCONSTANTS_H
