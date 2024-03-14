//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_BASEQUALITYFILTER_H
#define MUTECT2CPP_MASTER_BASEQUALITYFILTER_H

#include "HardFilter.h"


class BaseQualityFilter : public HardFilter{
private:
    double minMedianBaseQuality;

public:
    BaseQualityFilter(double minMedianBaseQuality);
    ErrorType errorType() override;
    bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine) override;
    std::vector<std::string> requiredAnnotations() override;
    std::string filterName() override;
    int filterIndex() override;

};


#endif //MUTECT2CPP_MASTER_BASEQUALITYFILTER_H
