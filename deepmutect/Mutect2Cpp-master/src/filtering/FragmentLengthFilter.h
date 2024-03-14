//
// Created by cluster on 22-11-15.
//

#ifndef MUTECT2CPP_MASTER_FRAGMENTLENGTHFILTER_H
#define MUTECT2CPP_MASTER_FRAGMENTLENGTHFILTER_H

#include "HardFilter.h"

class FragmentLengthFilter : public HardFilter{
private:
    double maxMedianFragmentLengthDifference;

public:
    explicit FragmentLengthFilter(double maxMedianFragmentLengthDifference);
    bool isArtifact(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine) override;
    ErrorType errorType() override;
    std::string filterName() override;
    std::vector<std::string> requiredAnnotations() override;
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_FRAGMENTLENGTHFILTER_H
