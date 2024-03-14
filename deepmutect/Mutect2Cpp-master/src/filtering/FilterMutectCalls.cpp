//
// Created by cluster on 22-11-10.
//

#include "FilterMutectCalls.h"

FilterMutectCalls::FilterMutectCalls(const std::string & normalSample) : MTFAC() ,filteringEngine(MTFAC, normalSample){

}

void FilterMutectCalls::nthPassApply(const std::shared_ptr<VariantContext> &vc,
                                     const std::shared_ptr<ReferenceContext> &referenceContext) {
    filteringEngine.accumulateData(vc, referenceContext);
}

void FilterMutectCalls::afterNthPass() {
    filteringEngine.learnParameters();
}

bool FilterMutectCalls::applyFilters(const std::shared_ptr<VariantContext> &vc,
                                     const std::shared_ptr<ReferenceContext> &referenceContext) {
    return filteringEngine.applyFiltersAndAccumulateOutputStats(vc, referenceContext);
}

