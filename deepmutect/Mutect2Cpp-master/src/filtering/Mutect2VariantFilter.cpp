//
// Created by cluster on 22-11-5.
//

#include "Mutect2VariantFilter.h"

double Mutect2VariantFilter::errorProbability(const std::shared_ptr<VariantContext> &vc,
                                              Mutect2FilteringEngine* filteringEngine,
                                              const std::shared_ptr<ReferenceContext> &referenceContext) {
    bool flag = true;
    for(const auto & str : requiredAnnotations()) {
        if(!vc->hasAttribute(str)) {
            flag = false;
            break;
        }
    }
    double result = flag ? calculateErrorProbability(vc, filteringEngine, referenceContext): 0;
    return Mutect2FilteringEngine::roundFinitePrecisionErrors(result);
}

void Mutect2VariantFilter::learnParameters() {

}

void Mutect2VariantFilter::learnParametersAndClearAccumulatedData() {
    learnParameters();
    clearAccumulatedData();
}

void Mutect2VariantFilter::clearAccumulatedData() {

}

void Mutect2VariantFilter::accumulateDataForLearning(const std::shared_ptr<VariantContext> &vc,
                                                     ErrorProbabilities errorProbabilities,
                                                     Mutect2FilteringEngine* filteringEngine) {

}

