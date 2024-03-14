//
// Created by cluster on 22-11-14.
//

#include "ReadPositionFilter.h"

ReadPositionFilter::ReadPositionFilter(double minMedianReadPosition) : minMedianReadPosition(minMedianReadPosition){

}

ErrorType ReadPositionFilter::errorType() {
    return ARTIFACT;
}

bool
ReadPositionFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    std::vector<int> readPositionByAllele;
    if(vc->hasAttribute(VCFConstants::MEDIAN_READ_POSITON_KEY)) {
        readPositionByAllele = vc->getAttributes().at(VCFConstants::MEDIAN_READ_POSITON_KEY).getAttributeAsIntVector();
    }
    return readPositionByAllele[0] > -1 && readPositionByAllele[0] < minMedianReadPosition;
}

std::string ReadPositionFilter::filterName() {
    return VCFConstants::READ_POSITION_FILTER_NAME;
}

std::vector<std::string> ReadPositionFilter::requiredAnnotations() {
    return {VCFConstants::MEDIAN_READ_POSITON_KEY};
}

int ReadPositionFilter::filterIndex() {
    return 10;
}

