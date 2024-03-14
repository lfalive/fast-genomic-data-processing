//
// Created by cluster on 22-11-10.
//

#ifndef MUTECT2CPP_MASTER_FILTERMUTECTCALLS_H
#define MUTECT2CPP_MASTER_FILTERMUTECTCALLS_H

#include "Mutect2FilteringEngine.h"


class FilterMutectCalls {
private:
    M2FiltersArgumentCollection MTFAC;
    Mutect2FilteringEngine filteringEngine;

public:
    FilterMutectCalls(const std::string & normalSample);
    void nthPassApply(const std::shared_ptr<VariantContext> & vc, const std::shared_ptr<ReferenceContext> & referenceContext);
    void afterNthPass();
    bool applyFilters(const std::shared_ptr<VariantContext> & vc, const std::shared_ptr<ReferenceContext> & referenceContext);
};


#endif //MUTECT2CPP_MASTER_FILTERMUTECTCALLS_H
