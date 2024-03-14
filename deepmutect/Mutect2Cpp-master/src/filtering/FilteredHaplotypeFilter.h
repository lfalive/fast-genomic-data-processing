//
// Created by cluster on 22-11-12.
//

#ifndef MUTECT2CPP_MASTER_FILTEREDHAPLOTYPEFILTER_H
#define MUTECT2CPP_MASTER_FILTEREDHAPLOTYPEFILTER_H

#include "Mutect2VariantFilter.h"

class FilteredHaplotypeFilter : public Mutect2VariantFilter{
private:
    std::optional<std::string> makePhasingString(std::shared_ptr<Genotype> &genotype);
    std::unordered_map<std::string, std::vector<std::pair<int, double>>> phasedProbabilities;
    std::unordered_map<std::string, std::vector<std::pair<int, double>>> accumulatingPhasedProbabilities;
    double maxIntraHaplotypeDistance;

public:
    FilteredHaplotypeFilter(double maxIntraHaplotypeDistance);
    virtual double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>);
    virtual void accumulateDataForLearning(const std::shared_ptr<VariantContext> & vc, ErrorProbabilities errorProbabilities, Mutect2FilteringEngine* filteringEngine);
    virtual void learnParameters();
    virtual void clearAccumulatedData();
    virtual std::string filterName();
    virtual ErrorType errorType();
    virtual std::vector<std::string> requiredAnnotations();
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_FILTEREDHAPLOTYPEFILTER_H
