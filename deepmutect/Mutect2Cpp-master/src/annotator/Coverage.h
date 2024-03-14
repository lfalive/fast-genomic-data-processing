//
// Created by lhh on 6/13/22.
//

#ifndef MUTECT2CPP_MASTER_COVERAGE_H
#define MUTECT2CPP_MASTER_COVERAGE_H

#include "PerAlleleAnnotation.h"

class Coverage : public InfoFieldAnnotation {
public:
    std::shared_ptr<std::map<std::string, AttributeValue>> annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

    std::vector<std::string> getKeyNames();
};


#endif //MUTECT2CPP_MASTER_COVERAGE_H
