//
// Created by 梦想家xixi on 2021/12/14.
//

#include "AssemblyRegionTrimmer.h"

AssemblyRegionTrimmer::AssemblyRegionTrimmer(ReadThreadingAssemblerArgumentCollection *assemblyArgs,
                                             SAMSequenceDictionary *sequenceDictionary, bool isGGA,
                                             bool emitReferenceConfidence) : assemblyArgs(assemblyArgs), sequenceDictionary(sequenceDictionary),
                                             emitReferenceConfidence(emitReferenceConfidence){
    Mutect2Utils::validateArg(assemblyArgs != nullptr, "null is not allowed there");
    checkUserArguments();
    usableExtension = isGGA ? assemblyArgs->ggaExtension : assemblyArgs->discoverExtension;
}

void AssemblyRegionTrimmer::checkUserArguments() {
    if(assemblyArgs->snpPadding < 0) {
        throw std::invalid_argument("paddingAroundSNPs");
    }
    if(assemblyArgs->indelPadding < 0) {
        throw std::invalid_argument("paddingAroundIndels");
    }
    if(assemblyArgs->discoverExtension < 0) {
        throw std::invalid_argument("maxDiscARExtension");
    }
    if(assemblyArgs->ggaExtension < 0) {
        throw std::invalid_argument("maxGGAAREExtension");
    }
}

std::shared_ptr<AssemblyRegionTrimmer_Result> AssemblyRegionTrimmer::trim(const std::shared_ptr<AssemblyRegion>& originalRegion,
                                                          std::set<std::shared_ptr<VariantContext>, VariantContextComparator> & allVariantsWithinExtendedRegion) {
    if(allVariantsWithinExtendedRegion.empty()) {
        return AssemblyRegionTrimmer_Result::noVariation(emitReferenceConfidence, originalRegion, assemblyArgs->snpPadding, usableExtension);
    }
    std::vector<std::shared_ptr<VariantContext>> withinActiveRegion;
	withinActiveRegion.reserve(allVariantsWithinExtendedRegion.size());
    std::shared_ptr<SimpleInterval> originalRegionRange = originalRegion->getSpan();
    bool foundNonSnp = false;
    std::shared_ptr<SimpleInterval> variantSpan = nullptr;
    for(const std::shared_ptr<VariantContext>& vc : allVariantsWithinExtendedRegion) {
        std::shared_ptr<SimpleInterval> vcLoc = std::make_shared<SimpleInterval>(vc->getContig(), vc->getStart(), vc->getEnd());
        if(originalRegionRange->overlaps(vcLoc)) {
            foundNonSnp = foundNonSnp || !vc->isSNP();
            variantSpan = variantSpan == nullptr ? vcLoc : variantSpan->spanWith(vcLoc);
            withinActiveRegion.emplace_back(vc);
        }
    }
    int padding = foundNonSnp ? assemblyArgs->indelPadding : assemblyArgs->snpPadding;
    if(variantSpan == nullptr) {
        return AssemblyRegionTrimmer_Result::noVariation(emitReferenceConfidence, originalRegion, padding, usableExtension);
    }

    if(assemblyArgs->dontTrimActiveRegions) {
        return AssemblyRegionTrimmer_Result::noTrimming(emitReferenceConfidence, originalRegion, padding, usableExtension, &withinActiveRegion);
    }
    std::shared_ptr<SimpleInterval> maximumSpan = originalRegionRange->expandWithinContig(usableExtension, sequenceDictionary);
    std::shared_ptr<SimpleInterval> idealSpan = variantSpan->expandWithinContig(padding, sequenceDictionary);
    std::shared_ptr<SimpleInterval> finalSpan = maximumSpan->intersect(idealSpan)->mergeWithContiguous(variantSpan);
    std::shared_ptr<SimpleInterval> callableSpan = emitReferenceConfidence ? variantSpan->intersect(originalRegionRange) : variantSpan;
    std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>  nonVariantRegions = nonVariantTargetRegions(originalRegion, callableSpan);

    return std::make_shared<AssemblyRegionTrimmer_Result>(emitReferenceConfidence, true, originalRegion, padding, usableExtension, &withinActiveRegion, nonVariantRegions, finalSpan, idealSpan, maximumSpan, variantSpan);
}

std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>
AssemblyRegionTrimmer::nonVariantTargetRegions(std::shared_ptr<AssemblyRegion> targetRegion, std::shared_ptr<SimpleInterval> variantSpan) {
    std::shared_ptr<SimpleInterval> targetRegionRange = targetRegion->getSpan();
    int finalStart = variantSpan->getStart();
    int finalStop = variantSpan->getEnd();
    int targetStart = targetRegionRange->getStart();
    int targetStop = targetRegionRange->getEnd();
    bool preTrimmingRequired = targetStart < finalStart;
    bool postTrimmingRequired = targetStop > finalStop;
    if(preTrimmingRequired) {
        std::string contig = targetRegionRange->getContig();
        return postTrimmingRequired ? std::make_shared<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>(std::make_shared<SimpleInterval>(contig, targetStart, finalStart-1), std::make_shared<SimpleInterval>(contig, finalStop+1, targetStop)) :
               std::make_shared<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>(std::make_shared<SimpleInterval>(contig, targetStart, finalStart-1), nullptr);
    } else if (postTrimmingRequired) {
        return std::make_shared<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>(nullptr, std::make_shared<SimpleInterval>(targetRegionRange->getContig(), finalStop+1, targetStop));
    } else {
        return std::make_shared<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>(nullptr, nullptr);
    }
}
