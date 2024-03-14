//
// Created by lhh on 5/20/22.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPEANNOTATION_H
#define MUTECT2CPP_MASTER_GENOTYPEANNOTATION_H

#include "VariantAnnotation.h"
#include "engine/ReferenceContext.h"
#include "VariantContext.h"
#include "Genotype.h"
#include "variantcontext/builder/GenotypeBuilder.h"
#include "utils/genotyper/AlleleLikelihoods.h"

/**
 * Represents an annotation that is computed for a single genotype.
 */
class GenotypeAnnotation : public VariantAnnotation{
public:
    const static std::string GENOTYPE_ANNOTATION;

    virtual void annotate(ReferenceContext& ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder& gb, AlleleLikelihoods<SAMRecord, Allele>* likelihoods) = 0;

    std::string toString();

    
};


#endif //MUTECT2CPP_MASTER_GENOTYPEANNOTATION_H
