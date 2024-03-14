//
// Created by cluster on 22-11-8.
//

#include "M2FiltersArgumentCollection.h"
#include "MathUtils.h"
#include <cfloat>

const double M2FiltersArgumentCollection:: DEFAULT_LOG_INDEL_PRIOR = MathUtils::log10ToLog(-7);
const double M2FiltersArgumentCollection:: DEFAULT_LOG_SNV_PRIOR = MathUtils::log10ToLog(-6);
const double M2FiltersArgumentCollection:: DEFAULT_LOG_INDEL_PRIOR_FOR_MITO = MathUtils::log10ToLog(-3.75);
const double M2FiltersArgumentCollection:: DEFAULT_LOG_SNV_PRIOR_FOR_MITO = MathUtils::log10ToLog(-2.5);
const double M2FiltersArgumentCollection:: DEFAULT_INITIAL_POSTERIOR_THRESHOLD = 0.1;
const double M2FiltersArgumentCollection:: DEFAULT_MAX_FALSE_DISCOVERY_RATE = 0.05;
const double M2FiltersArgumentCollection:: DEFAULT_F_SCORE_BETA = 1.0;
const double M2FiltersArgumentCollection:: DEFAULT_INITIAL_LOG_PRIOR_OF_VARIANT_VERSUS_ARTIFACT = MathUtils::log10ToLog(-1);
const double M2FiltersArgumentCollection::DEFAULT_MAX_N_RATIO = 1.0 / 0.0;

double M2FiltersArgumentCollection::getLogIndelPrior() {
    return mitochondria && logIndelPrior == DEFAULT_LOG_INDEL_PRIOR ? DEFAULT_LOG_INDEL_PRIOR_FOR_MITO : logIndelPrior;
}

M2FiltersArgumentCollection::M2FiltersArgumentCollection() {
    mitochondria = false;
    logIndelPrior = DEFAULT_LOG_INDEL_PRIOR;
    logSNVPrior = DEFAULT_LOG_SNV_PRIOR;
    fScoreBeta = DEFAULT_F_SCORE_BETA;
    maxFalsePositiveRate = DEFAULT_MAX_FALSE_DISCOVERY_RATE;
    initialPosteriorThreshold = DEFAULT_INITIAL_POSTERIOR_THRESHOLD;
    initialLogPriorOfVariantVersusArtifact = DEFAULT_INITIAL_LOG_PRIOR_OF_VARIANT_VERSUS_ARTIFACT;
    minMedianBaseQuality = DEFAULT_MIN_MEDIAN_BASE_QUALITY;
    nRatio = DEFAULT_MAX_N_RATIO;
    minMedianReadPosition = DEFAULT_MIN_MEDIAN_READ_POSITION;
}

double M2FiltersArgumentCollection::getLogSnvPrior() {
    return mitochondria && logSNVPrior == DEFAULT_LOG_SNV_PRIOR ? DEFAULT_LOG_SNV_PRIOR_FOR_MITO : logSNVPrior;
}
