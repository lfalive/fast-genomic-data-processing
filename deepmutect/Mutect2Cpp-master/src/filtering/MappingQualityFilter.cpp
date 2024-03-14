//
// Created by cluster on 22-11-14.
//

#include "MappingQualityFilter.h"

MappingQualityFilter::MappingQualityFilter(double minMedianMappingQuality, int longIndelSize) : minMedianMappingQuality(minMedianMappingQuality), longIndelSize(longIndelSize){

}

ErrorType MappingQualityFilter::errorType() {
    return ARTIFACT;
}

bool
MappingQualityFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    std::vector<int> indelLengths = vc->getIndelLengths();
    int indelLength = 0;
    if(!indelLengths.empty()) {
        for(int i : indelLengths) {
            if(indelLength < abs(i)) {
                indelLength = abs(i);
            }
        }
    }
    std::vector<int> mappingQualityByAllele = vc->getAttributeAsIntVector(VCFConstants::MEDIAN_MAPPING_QUALITY_KEY, {});
    return mappingQualityByAllele[indelLength < longIndelSize ? 1 : 0] < minMedianMappingQuality;
}

std::string MappingQualityFilter::filterName() {
    return VCFConstants::MEDIAN_MAPPING_QUALITY_FILTER_NAME;
}

std::vector<std::string> MappingQualityFilter::requiredAnnotations() {
    return {VCFConstants::MEDIAN_MAPPING_QUALITY_KEY};
}

int MappingQualityFilter::filterIndex() {
    return 19;
}
