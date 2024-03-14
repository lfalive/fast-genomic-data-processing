//
// Created by cluster on 22-11-14.
//

#include "DuplicatedAltReadFilter.h"

DuplicatedAltReadFilter::DuplicatedAltReadFilter(int uniqueAltReadCount) : uniqueAltReadCount(uniqueAltReadCount){

}

ErrorType DuplicatedAltReadFilter::errorType() {
    return ARTIFACT;
}

bool DuplicatedAltReadFilter::isArtifact(const std::shared_ptr<VariantContext> &vc,
                                         Mutect2FilteringEngine *filteringEngine) {
    return vc->getAttributeAsInt(VCFConstants::UNIQUE_ALT_READ_SET_COUNT_KEY, 1) <= uniqueAltReadCount;
}

std::string DuplicatedAltReadFilter::filterName() {
    return VCFConstants::DUPLICATED_EVIDENCE_FILTER_NAME;
}

std::vector<std::string> DuplicatedAltReadFilter::requiredAnnotations() {
    return {VCFConstants::UNIQUE_ALT_READ_SET_COUNT_KEY};
}

int DuplicatedAltReadFilter::filterIndex() {
    return 18;
}
