//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_NORMALARTIFACTFILTER_H
#define MUTECT2CPP_MASTER_NORMALARTIFACTFILTER_H

#include "Mutect2VariantFilter.h"


class NormalArtifactFilter : public Mutect2VariantFilter{
private:
    constexpr static double MIN_NORMAL_ARTIFACT_RATIO = 0.1;
    const static int IMPUTED_NORMAL_BASE_QUALITY = 30;

public:
    ErrorType errorType() override;
    double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>) override;
    std::string filterName() override;
    std::vector<std::string> requiredAnnotations() override;
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_NORMALARTIFACTFILTER_H
