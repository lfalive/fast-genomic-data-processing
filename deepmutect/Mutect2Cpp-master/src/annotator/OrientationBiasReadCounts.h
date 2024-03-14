//
// Created by lhh on 6/7/22.
//

#ifndef MUTECT2CPP_MASTER_ORIENTATIONBIASREADCOUNTS_H
#define MUTECT2CPP_MASTER_ORIENTATIONBIASREADCOUNTS_H

#include "GenotypeAnnotation.h"

class OrientationBiasReadCounts : public GenotypeAnnotation {
private:
    const static int MINIMUM_BASE_QUALITY = 20;

    static bool isUsableRead(shared_ptr<SAMRecord> read);

    static int getReadBaseQuality(shared_ptr<SAMRecord> read, int refLoc);

public:
    void annotate(ReferenceContext &ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder &gb,
             AlleleLikelihoods<SAMRecord, Allele> *likelihoods);

    std::vector<std::string> getKeyNames();
};


#endif //MUTECT2CPP_MASTER_ORIENTATIONBIASREADCOUNTS_H
