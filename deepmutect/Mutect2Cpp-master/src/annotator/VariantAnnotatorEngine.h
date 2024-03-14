//
// Created by lhh on 5/20/22.
//

#ifndef MUTECT2CPP_MASTER_VARIANTANNOTATORENGINE_H
#define MUTECT2CPP_MASTER_VARIANTANNOTATORENGINE_H

#include <vector>
#include <set>
#include <string>
#include <memory>
#include "VariantContext.h"
#include "utils/genotyper/AlleleLikelihoods.h"
#include "engine/ReferenceContext.h"
#include "InfoFieldAnnotation.h"
#include "GenotypeAnnotation.h"

/**
 * The class responsible for computing annotations for variants.
 * Annotations are auto-discovered - ie, any class that extends {@link VariantAnnotation} and
 * lives in this package is treated as an annotation and the engine will attempt to create instances of it
 * by calling the non-arg constructor (loading will fail if there is no no-arg constructor).
 */
class VariantAnnotatorEngine {
private:
    std::vector<shared_ptr<InfoFieldAnnotation>> infoAnnotations;
    std::vector<shared_ptr<GenotypeAnnotation>> genotypeAnnotations;
    std::set<std::string> reducibleKeys;
    bool useRawAnnotations;
    bool keepRawCombinedAnnotations;

public:
    VariantAnnotatorEngine();

    VariantAnnotatorEngine(std::vector<shared_ptr<InfoFieldAnnotation>>&& InfoFieldAnnotationList, std::vector<shared_ptr<GenotypeAnnotation>>&& GenotypeAnnotationList);

    shared_ptr<VariantContext> annotateContext(shared_ptr<VariantContext> vc, ReferenceContext& ref, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

    shared_ptr<GenoTypesContext> annotateGenotypes(ReferenceContext& ref, shared_ptr<VariantContext> vc,  AlleleLikelihoods<SAMRecord, Allele>* likelihoods);
};


#endif //MUTECT2CPP_MASTER_VARIANTANNOTATORENGINE_H
