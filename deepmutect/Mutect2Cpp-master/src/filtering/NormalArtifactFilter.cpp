//
// Created by cluster on 22-11-14.
//

#include "NormalArtifactFilter.h"
#include "MathUtils.h"
#include "clustering/BinomialDistribution.h"

ErrorType NormalArtifactFilter::errorType() {
    return ARTIFACT;
}

double NormalArtifactFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                       Mutect2FilteringEngine *filteringEngine,
                                                       std::shared_ptr<ReferenceContext>) {
    std::vector<double> tumorLods = Mutect2FilteringEngine::getTumorLogOdds(vc);
    int indexOfMaxTumorLod = MathUtils::maxElementIndex(tumorLods);

    std::vector<int> tumorAlleleDepths = filteringEngine->sumADsOverSamples(vc, true, false);
    int tumorDepth = 0;
    for(auto i : tumorAlleleDepths) {
        tumorDepth += i;
    }
    int tumorAltDepth = tumorAlleleDepths[indexOfMaxTumorLod + 1];

    std::vector<int> normalAlleleDepths = filteringEngine->sumADsOverSamples(vc, false, true);
    int normalDepth = 0;
    for(auto i : normalAlleleDepths) {
        normalDepth += i;
    }
    int normalAltDepth = normalAlleleDepths[indexOfMaxTumorLod + 1];

    double tumorAlleleFraction = (double) tumorAltDepth / tumorDepth;
    double normalAlleleFraction = normalDepth == 0 ? 0 : (double) normalAltDepth / normalDepth;

    if (normalAlleleFraction < MIN_NORMAL_ARTIFACT_RATIO * tumorAlleleFraction)  {
        return 0.0;
    }

    std::shared_ptr<std::vector<double>> normalArtifactNegativeLogOdds = MathUtils::applyToArrayInPlace(std::make_shared<std::vector<double>>(vc->getAttributes().at(VCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY).getAttributeAsDoubleVector()), MathUtils::log10ToLog);
    double normalArtifactProbability = filteringEngine->posteriorProbabilityOfNormalArtifact((*normalArtifactNegativeLogOdds)[indexOfMaxTumorLod]);

    int medianRefBaseQuality = IMPUTED_NORMAL_BASE_QUALITY;
    if(vc->hasAttribute(VCFConstants::MEDIAN_BASE_QUALITY_KEY)) {
        medianRefBaseQuality = vc->getAttributes().at(VCFConstants::MEDIAN_BASE_QUALITY_KEY).getAttributeAsIntVector()[0];
    }
    double normalPValue = 1 - BinomialDistribution(normalDepth, QualityUtils::qualToErrorProb((double)medianRefBaseQuality)).cumulativeProbability(normalAltDepth - 1);

    return normalPValue < M2FiltersArgumentCollection::normalPileupPValueThreshold ? 1.0 : normalArtifactProbability;
}

std::string NormalArtifactFilter::filterName() {
    return VCFConstants::ARTIFACT_IN_NORMAL_FILTER_NAME;
}

std::vector<std::string> NormalArtifactFilter::requiredAnnotations() {
    return {VCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY, VCFConstants::TUMOR_LOG_10_ODDS_KEY};
}

int NormalArtifactFilter::filterIndex() {
    return 21;
}
