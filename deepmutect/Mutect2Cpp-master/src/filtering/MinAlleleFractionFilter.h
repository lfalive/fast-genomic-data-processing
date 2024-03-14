//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_MINALLELEFRACTIONFILTER_H
#define MUTECT2CPP_MASTER_MINALLELEFRACTIONFILTER_H

#include "HardFilter.h"

class MinAlleleFractionFilter : public HardFilter{
private:
    double minAf;

public:
    explicit MinAlleleFractionFilter(double minAf);
    ErrorType errorType() override;
    bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine) override;
    std::string filterName() override;
    std::vector<std::string> requiredAnnotations() override;
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_MINALLELEFRACTIONFILTER_H
