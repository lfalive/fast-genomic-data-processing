//
// Created by cluster on 22-11-14.
//

#include "HardFilter.h"

double HardFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                             Mutect2FilteringEngine *filteringEngine,
                                             std::shared_ptr<ReferenceContext>) {
    return isArtifact(vc, filteringEngine) ? 1 : 0;
}
