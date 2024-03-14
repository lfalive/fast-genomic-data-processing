//
// Created by cluster on 22-11-14.
//

#include "MinAlleleFractionFilter.h"

MinAlleleFractionFilter::MinAlleleFractionFilter(double minAf) : minAf(minAf){

}

ErrorType MinAlleleFractionFilter::errorType() {
    return ARTIFACT;
}

bool MinAlleleFractionFilter::isArtifact(const std::shared_ptr<VariantContext> &vc,
                                         Mutect2FilteringEngine *filteringEngine) {
    for(auto & k : *vc->getGenotypes()->getGenotypes()) {
        if(filteringEngine->isTumor(k.get())) {
            if(k->hasExtendedAttribute(VCFConstants::ALLELE_FREQUENCY_KEY)) {
                std::vector<double> alleleFractions;
                if(k->hasExtendedAttribute(VCFConstants::ALLELE_FRACTION_KEY)) {
                    alleleFractions = k->getExtendedAttribute(
                            VCFConstants::ALLELE_FRACTION_KEY).getAttributeAsDoubleVector();
                } else {
                    alleleFractions = {1};
                }
                int numRealAlleles = vc->hasSymbolicAlleles() ? alleleFractions.size() - 1 : alleleFractions.size();
                double max = alleleFractions[0];
                for(auto i : alleleFractions) {
                    if(max < i) {
                        max = i;
                    }
                }
                return max < minAf;
            }
        }
    }
    return false;
}

std::string MinAlleleFractionFilter::filterName() {
    return VCFConstants::ALLELE_FRACTION_FILTER_NAME;
}

std::vector<std::string> MinAlleleFractionFilter::requiredAnnotations() {
    return {};
}

int MinAlleleFractionFilter::filterIndex() {
    return 1;
}
