//
// Created by lhh on 6/7/22.
//

#ifndef MUTECT2CPP_MASTER_DEPTHPERSAMPLEHC_H
#define MUTECT2CPP_MASTER_DEPTHPERSAMPLEHC_H

#include "GenotypeAnnotation.h"

class DepthPerSampleHC : public GenotypeAnnotation{
public:
    std::vector<std::string> getKeyNames();

    void annotate(ReferenceContext& ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder& gb, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

};


#endif //MUTECT2CPP_MASTER_DEPTHPERSAMPLEHC_H
