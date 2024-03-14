//
// Created by cluster on 22-11-5.
//

#ifndef MUTECT2CPP_MASTER_TUMOREVIDENCEFILTER_H
#define MUTECT2CPP_MASTER_TUMOREVIDENCEFILTER_H

#include "Mutect2VariantFilter.h"

class TumorEvidenceFilter : public Mutect2VariantFilter{
public:
    virtual double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>);
    virtual std::vector<std::string> requiredAnnotations();
    virtual ErrorType errorType();
    virtual std::string filterName();
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_TUMOREVIDENCEFILTER_H
