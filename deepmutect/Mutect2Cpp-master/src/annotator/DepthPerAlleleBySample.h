//
// Created by lhh on 6/6/22.
//

#ifndef MUTECT2CPP_MASTER_DEPTHPERALLELEBYSAMPLE_H
#define MUTECT2CPP_MASTER_DEPTHPERALLELEBYSAMPLE_H

#include "GenotypeAnnotation.h"

class DepthPerAlleleBySample : public GenotypeAnnotation {
public:
    void annotate(ReferenceContext& ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder& gb, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

	std::vector<int> annotateWithLikelihoods(shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, vector<shared_ptr<Allele>>& alleles, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

    std::vector<std::string> getKeyNames();
};


#endif //MUTECT2CPP_MASTER_DEPTHPERALLELEBYSAMPLE_H
