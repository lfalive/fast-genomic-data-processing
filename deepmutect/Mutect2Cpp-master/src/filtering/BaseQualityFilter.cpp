//
// Created by cluster on 22-11-14.
//

#include "BaseQualityFilter.h"
#include "Mutect2FilteringEngine.h"
#include "MathUtils.h"

BaseQualityFilter::BaseQualityFilter(double minMedianBaseQuality) : minMedianBaseQuality(minMedianBaseQuality){

}

ErrorType BaseQualityFilter::errorType() {
    return ARTIFACT;
}

bool BaseQualityFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    std::vector<int> baseQualityByAllele = vc->getAttributeAsIntVector(VCFConstants::MEDIAN_BASE_QUALITY_KEY, {});
    std::vector<double> tumorLods = Mutect2FilteringEngine::getTumorLogOdds(vc);
    int indexOfMaxTumorLod = MathUtils::maxElementIndex(tumorLods);

    return baseQualityByAllele[indexOfMaxTumorLod + 1] < minMedianBaseQuality;
}

std::vector<std::string> BaseQualityFilter::requiredAnnotations() {
    return {VCFConstants::MEDIAN_BASE_QUALITY_KEY};
}

std::string BaseQualityFilter::filterName() {
    return VCFConstants::MEDIAN_BASE_QUALITY_FILTER_NAME;
}

int BaseQualityFilter::filterIndex() {
    return 14;
}
