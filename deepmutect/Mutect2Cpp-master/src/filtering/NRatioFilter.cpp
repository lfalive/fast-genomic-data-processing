//
// Created by cluster on 22-11-14.
//

#include "NRatioFilter.h"

NRatioFilter::NRatioFilter(double maxRatio) : maxNRatio(maxRatio){

}

ErrorType NRatioFilter::errorType() {
    return ARTIFACT;
}

bool NRatioFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    std::vector<int> ADs = filteringEngine->sumADsOverSamples(vc, true, true);
    int altCount = 0;
    for(auto i : ADs) {
        altCount += i;
    }
    altCount -= ADs[0];
    if (altCount == 0 ) {
        return false;
    }

    int NCount = vc->getAttributeAsInt(VCFConstants::N_COUNT_KEY, 0);

    return (double) NCount / altCount >= maxNRatio;
}

std::string NRatioFilter::filterName() {
    return VCFConstants::N_RATIO_FILTER_NAME;
}

std::vector<std::string> NRatioFilter::requiredAnnotations() {
    return {VCFConstants::N_COUNT_KEY};
}

int NRatioFilter::filterIndex() {
    return 5;
}
