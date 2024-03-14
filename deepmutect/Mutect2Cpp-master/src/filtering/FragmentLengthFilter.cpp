//
// Created by cluster on 22-11-15.
//

#include "FragmentLengthFilter.h"

FragmentLengthFilter::FragmentLengthFilter(double maxMedianFragmentLengthDifference) : maxMedianFragmentLengthDifference(maxMedianFragmentLengthDifference){

}

bool
FragmentLengthFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    std::vector<int> fragmentLengthByAllele = {0};
    if(vc->hasAttribute(VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY)) {
        fragmentLengthByAllele = vc->getAttributes().at(VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY).getAttributeAsIntVector();
    }
    return std::abs(fragmentLengthByAllele[1] - fragmentLengthByAllele[0]) > maxMedianFragmentLengthDifference;
}

ErrorType FragmentLengthFilter::errorType() {
    return ARTIFACT;
}

std::string FragmentLengthFilter::filterName() {
    return {VCFConstants::MEDIAN_FRAGMENT_LENGTH_DIFFERENCE_FILTER_NAME};
}

std::vector<std::string> FragmentLengthFilter::requiredAnnotations() {
    return {VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY};
}

int FragmentLengthFilter::filterIndex() {
    return 13;
}
