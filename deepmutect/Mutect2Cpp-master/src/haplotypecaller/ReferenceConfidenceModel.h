//
// Created by 梦想家xixi on 2022/1/15.
//

#ifndef MUTECT2CPP_MASTER_REFERENCECONFIDENCEMODEL_H
#define MUTECT2CPP_MASTER_REFERENCECONFIDENCEMODEL_H

#include "Haplotype.h"
#include "AssemblyRegion.h"

class ReferenceConfidenceModel {
public:
    static std::shared_ptr<Haplotype> createReferenceHaplotype(const std::shared_ptr<AssemblyRegion>& activeRegion, const std::shared_ptr<uint8_t[]>& refBase, int &length,const std::shared_ptr<SimpleInterval> & paddedReferenceLoc);
};


#endif //MUTECT2CPP_MASTER_REFERENCECONFIDENCEMODEL_H
