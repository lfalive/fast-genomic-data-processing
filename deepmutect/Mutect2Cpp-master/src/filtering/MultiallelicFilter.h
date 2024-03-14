//
// Created by cluster on 22-11-15.
//

#ifndef MUTECT2CPP_MASTER_MULTIALLELICFILTER_H
#define MUTECT2CPP_MASTER_MULTIALLELICFILTER_H

#include "HardFilter.h"

class MultiallelicFilter : public HardFilter{
private:
    int numAltAllelesThreshold;
    constexpr static double MULTIALLELIC_LOD_THRESHOLD = 5.0;

public:
    explicit MultiallelicFilter(int numAltAllelesThreshold);
    ErrorType errorType() override;
    bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine) override;
    std::vector<std::string> requiredAnnotations() override;
    std::string filterName() override;
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_MULTIALLELICFILTER_H
