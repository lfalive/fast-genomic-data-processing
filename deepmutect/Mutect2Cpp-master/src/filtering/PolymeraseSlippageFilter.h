//
// Created by cluster on 22-11-15.
//

#ifndef MUTECT2CPP_MASTER_POLYMERASESLIPPAGEFILTER_H
#define MUTECT2CPP_MASTER_POLYMERASESLIPPAGEFILTER_H

#include "Mutect2VariantFilter.h"

class PolymeraseSlippageFilter : public Mutect2VariantFilter{
private:
    int minSlippageLength;
    double slippageRate;

public:
    PolymeraseSlippageFilter(int minSlippageLength, double slippageRate);
    ErrorType errorType() override;
    double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>) override;
    std::string filterName() override;
    std::vector<std::string> requiredAnnotations() override;
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_POLYMERASESLIPPAGEFILTER_H
