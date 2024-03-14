//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_DUPLICATEDALTREADFILTER_H
#define MUTECT2CPP_MASTER_DUPLICATEDALTREADFILTER_H

#include "HardFilter.h"

class DuplicatedAltReadFilter : public HardFilter{
private:
    int uniqueAltReadCount;

public:
    explicit DuplicatedAltReadFilter(int uniqueAltReadCount);
    ErrorType errorType() override;
    bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine) override;
    std::string filterName() override;
    std::vector<std::string> requiredAnnotations() override;
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_DUPLICATEDALTREADFILTER_H
