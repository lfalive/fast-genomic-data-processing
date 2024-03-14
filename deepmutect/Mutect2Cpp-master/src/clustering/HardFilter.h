//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_HARDFILTER_H
#define MUTECT2CPP_MASTER_HARDFILTER_H

#include "Mutect2VariantFilter.h"

class HardFilter : public Mutect2VariantFilter{
public:
    virtual double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>);
    virtual bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine) = 0;
};


#endif //MUTECT2CPP_MASTER_HARDFILTER_H
