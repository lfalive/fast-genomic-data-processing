//
// Created by cluster on 22-11-14.
//

#include "PanelOfNormalsFilter.h"

ErrorType PanelOfNormalsFilter::errorType() {
    return ARTIFACT;
}

bool
PanelOfNormalsFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    return vc->hasAttribute(VCFConstants::IN_PON_KEY);
}

std::string PanelOfNormalsFilter::filterName() {
    return VCFConstants::PON_FILTER_NAME;
}

std::vector<std::string> PanelOfNormalsFilter::requiredAnnotations() {
    return {};
}

int PanelOfNormalsFilter::filterIndex() {
    return 2;
}
