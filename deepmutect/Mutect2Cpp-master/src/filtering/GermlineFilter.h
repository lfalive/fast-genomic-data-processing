//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_GERMLINEFILTER_H
#define MUTECT2CPP_MASTER_GERMLINEFILTER_H

#include "Mutect2VariantFilter.h"

class GermlineFilter : public Mutect2VariantFilter{
private:
    constexpr static double MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT = 0.9;
    constexpr static double EPSILON = 1.0e-10;

    double computeMinorAlleleFraction(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, const std::vector<int> & alleleCounts);

public:
    double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>) override;

    static double germlineProbability(double normalLogOdds, double logOddsOfGermlineHetVsSomatic, double logOddsOfGermlineHomAltVsSomatic, double populationAF, double logPriorSomatic);

    std::string filterName() override;

    std::vector<std::string> requiredAnnotations() override;

    ErrorType errorType() override;

    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_GERMLINEFILTER_H
