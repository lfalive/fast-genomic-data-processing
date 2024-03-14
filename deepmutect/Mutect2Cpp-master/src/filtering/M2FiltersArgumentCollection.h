//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H


class
M2FiltersArgumentCollection {
public:
    double getLogIndelPrior();
    double getLogSnvPrior();
    bool mitochondria;
    double logIndelPrior;
    double logSNVPrior;
    double initialPosteriorThreshold;
    double fScoreBeta;
    double maxFalsePositiveRate;
    double initialLogPriorOfVariantVersusArtifact;
    double nRatio;
    double minAf = 0;
    int minMedianBaseQuality;
    int minMedianReadPosition;
    int minMedianMappingQuality = DEFAULT_MIN_MEDIAN_MAPPING_QUALITY;
    int longIndelLength = DEFAULT_LONG_INDEL_SIZE;
    int uniqueAltReadCount = DEFAULT_MIN_UNIQUE_ALT_READS;
    int minSlippageLength = DEFAULT_MIN_SLIPPAGE_LENGTH;
    double slippageRate = DEFAULT_SLIPPAGE_RATE;
    constexpr static double normalPileupPValueThreshold = 0.0001;

    M2FiltersArgumentCollection();

private:
    static const double  DEFAULT_LOG_INDEL_PRIOR;
    static const double  DEFAULT_LOG_SNV_PRIOR;
    static const double  DEFAULT_LOG_INDEL_PRIOR_FOR_MITO;
    static const double  DEFAULT_LOG_SNV_PRIOR_FOR_MITO;
    static const double  DEFAULT_INITIAL_POSTERIOR_THRESHOLD;
    static const double  DEFAULT_F_SCORE_BETA;
    static const double  DEFAULT_MAX_FALSE_DISCOVERY_RATE;
    static const double  DEFAULT_INITIAL_LOG_PRIOR_OF_VARIANT_VERSUS_ARTIFACT;
    static const int DEFAULT_MIN_MEDIAN_BASE_QUALITY = 20;
    static const int DEFAULT_MIN_MEDIAN_MAPPING_QUALITY = 30;
    static const int DEFAULT_LONG_INDEL_SIZE = 5;
    static const int DEFAULT_MIN_UNIQUE_ALT_READS = 0;
    static const double DEFAULT_MAX_N_RATIO;
    constexpr static double DEFAULT_SLIPPAGE_RATE = 0.1;
    static const int DEFAULT_MIN_MEDIAN_READ_POSITION = 1;
    static const int DEFAULT_MIN_SLIPPAGE_LENGTH = 8;
};


#endif //MUTECT2CPP_MASTER_M2FILTERSARGUMENTCOLLECTION_H
