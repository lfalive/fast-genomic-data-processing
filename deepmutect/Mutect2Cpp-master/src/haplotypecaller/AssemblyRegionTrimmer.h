//
// Created by 梦想家xixi on 2021/12/14.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_H
#define MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_H

#include "ReadThreadingAssemblerArgumentCollection.h"
#include "samtools/SAMSequenceDictionary.h"
#include "AssemblyRegionTrimmer_Result.h"

class AssemblyRegionTrimmer {
private:
    int usableExtension;
    bool emitReferenceConfidence;
    ReadThreadingAssemblerArgumentCollection* assemblyArgs;
    SAMSequenceDictionary* sequenceDictionary;
    void checkUserArguments();

public:
    AssemblyRegionTrimmer(ReadThreadingAssemblerArgumentCollection* assemblyArgs, SAMSequenceDictionary* sequenceDictionary, bool isGGA, bool emitReferenceConfidence);
    std::shared_ptr<AssemblyRegionTrimmer_Result> trim(const std::shared_ptr<AssemblyRegion>& originalRegion, std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & allVariantsWithinExtendedRegion);
    std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>> nonVariantTargetRegions(std::shared_ptr<AssemblyRegion> targetRegion, std::shared_ptr<SimpleInterval> variantSpan);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_H
