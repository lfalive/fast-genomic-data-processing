//
// Created by cluster on 22-11-5.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2VARIANTFILTER_H
#define MUTECT2CPP_MASTER_MUTECT2VARIANTFILTER_H

#include "variantcontext/VariantContext.h"
#include "Mutect2FilteringEngine.h"
#include "engine/ReferenceContext.h"
#include "ErrorType.h"

class ErrorProbabilities;
class Mutect2FilteringEngine;

class Mutect2VariantFilter {
public:
    Mutect2VariantFilter() = default;
    double errorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, const std::shared_ptr<ReferenceContext>& referenceContext);
    virtual void accumulateDataForLearning(const std::shared_ptr<VariantContext> & vc, ErrorProbabilities errorProbabilities, Mutect2FilteringEngine* filteringEngine);
    virtual ErrorType errorType() = 0;
    virtual void learnParametersAndClearAccumulatedData();
    virtual ~Mutect2VariantFilter() {}
    virtual std::string filterName() = 0;

    virtual int filterIndex() = 0;

protected:
    virtual std::vector<std::string> requiredAnnotations() = 0;

    virtual double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>) = 0;
    virtual void learnParameters();
    virtual void clearAccumulatedData();
};


#endif //MUTECT2CPP_MASTER_MUTECT2VARIANTFILTER_H
