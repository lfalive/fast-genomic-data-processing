//
// Created by 梦想家xixi on 2021/12/14.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H
#define MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H


#include "AssemblyRegion.h"
#include "VariantContext.h"

class AssemblyRegionTrimmer_Result {
protected:
    bool needsTrimming;
    std::shared_ptr<AssemblyRegion> originalRegion;
    std::shared_ptr<SimpleInterval> callableSpan;
    std::shared_ptr<SimpleInterval> maximumSpan;
    std::shared_ptr<SimpleInterval> extendedSpan;
    std::shared_ptr<SimpleInterval> idealSpan;
    std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>> nonVariantFlanks;
    std::vector<std::shared_ptr<VariantContext>> callableEvents;
    int padding;
    int usableExtension;
    std::shared_ptr<AssemblyRegion> callableRegion;

private:
    std::shared_ptr<AssemblyRegion> leftFlankRegion;
    std::shared_ptr<AssemblyRegion> rightFlankRegion;
    bool emitReferenceConfidence;

public:
    AssemblyRegionTrimmer_Result(bool emitReferenceConfidence, bool needsTrimming, std::shared_ptr<AssemblyRegion>  originalRegion, int padding, int extension,
                                 std::vector<std::shared_ptr<VariantContext>> * overlappingEvents, std::shared_ptr<std::pair<std::shared_ptr<SimpleInterval>, std::shared_ptr<SimpleInterval>>>  nonVariantFlanks,
                                 const std::shared_ptr<SimpleInterval>& extendedSpan, std::shared_ptr<SimpleInterval>  idealSpan, std::shared_ptr<SimpleInterval>  maximumSpan, const std::shared_ptr<SimpleInterval>& callableSpan);
    static std::shared_ptr<AssemblyRegionTrimmer_Result> noVariation(bool emitReferenceConfidence, const std::shared_ptr<AssemblyRegion>& targetRegion, int padding, int usableExtension);
    ~AssemblyRegionTrimmer_Result();
    static std::shared_ptr<AssemblyRegionTrimmer_Result> noTrimming(bool emitReferenceConfidence, const std::shared_ptr<AssemblyRegion>& targetRegion, int padding, int usableExtension, std::vector<std::shared_ptr<VariantContext>> * events);
    bool isVariationPresent();
    std::shared_ptr<AssemblyRegion> getCallableRegion();
    bool getNeedsTrimming() const;
	void printInfo();
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGIONTRIMMER_RESULT_H
