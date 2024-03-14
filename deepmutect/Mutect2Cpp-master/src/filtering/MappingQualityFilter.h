//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_MAPPINGQUALITYFILTER_H
#define MUTECT2CPP_MASTER_MAPPINGQUALITYFILTER_H

#include "HardFilter.h"

class MappingQualityFilter : public HardFilter{
private:
    double minMedianMappingQuality;
    int longIndelSize;

public:
    MappingQualityFilter(double minMedianMappingQuality, int longIndelSize);
    virtual std::string filterName();
    virtual std::vector<std::string> requiredAnnotations();
    ErrorType errorType() override;
    virtual bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine);
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_MAPPINGQUALITYFILTER_H
