//
// Created by cluster on 22-11-15.
//

#include "MultiallelicFilter.h"

MultiallelicFilter::MultiallelicFilter(int numAltAllelesThreshold) : numAltAllelesThreshold(numAltAllelesThreshold){

}

ErrorType MultiallelicFilter::errorType() {
    return ARTIFACT;
}

bool
MultiallelicFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    std::vector<double> tumorLods = Mutect2FilteringEngine::getTumorLogOdds(vc);
    long numPassingAltAlleles = 0;
    for(double i : tumorLods) {
        if(i > MULTIALLELIC_LOD_THRESHOLD) {
            numPassingAltAlleles++;
        }
    }
    return numPassingAltAlleles > numAltAllelesThreshold;
}

std::vector<std::string> MultiallelicFilter::requiredAnnotations() {
    return {VCFConstants::TUMOR_LOG_10_ODDS_KEY};
}

std::string MultiallelicFilter::filterName() {
    return VCFConstants::MULTIALLELIC_FILTER_NAME;
}

int MultiallelicFilter::filterIndex() {
    return 12;
}
