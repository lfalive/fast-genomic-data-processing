//
// Created by cluster on 22-11-9.
//

#include "ErrorProbabilities.h"

ErrorProbabilities::ErrorProbabilities(std::vector<Mutect2VariantFilter *> filters,
                                       const std::shared_ptr<VariantContext> &vc,
                                       Mutect2FilteringEngine* filteringEngine, const std::shared_ptr<ReferenceContext>& referenceContext) {
    for(auto filter : filters){
        probabilitiesByFilter.insert({filter, filter->errorProbability(vc, filteringEngine, referenceContext)});
    }
    probabilitiesByType.insert({ARTIFACT, 0.0});
    probabilitiesByType.insert({NON_SOMATIC, 0.0});
    probabilitiesByType.insert({SEQUENCING, 0.0});
    for(auto kv : probabilitiesByFilter) {
        if(kv.second > probabilitiesByType.at(kv.first->errorType())) {
            probabilitiesByType[kv.first->errorType()] = kv.second;
        }
    }
    double trueProbability = 1;
    for(auto kv : probabilitiesByType) {
        trueProbability *= (1 - kv.second);
    }

    errorProbability = Mutect2FilteringEngine::roundFinitePrecisionErrors(1 - trueProbability);
}
