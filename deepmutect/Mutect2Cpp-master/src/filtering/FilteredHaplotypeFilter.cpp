//
// Created by cluster on 22-11-12.
//

#include "FilteredHaplotypeFilter.h"

double FilteredHaplotypeFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                          Mutect2FilteringEngine *filteringEngine,
                                                          std::shared_ptr<ReferenceContext>) {
    std::vector<std::shared_ptr<Genotype>> genotypes;
    for(auto & k : *vc->getGenotypes()->getGenotypes()) {
        if(filteringEngine->isTumor(k.get())) {
            genotypes.emplace_back(k);
        }
    }
    double max = 0.0;
    std::shared_ptr<Genotype> tumorGenotype = genotypes[0];
    for(auto & k : genotypes) {
        if(k->hasExtendedAttribute(VCFConstants::ALLELE_FREQUENCY_KEY)) {
            std::vector<double> tmp = tumorGenotype->getExtendedAttribute(VCFConstants::ALLELE_FREQUENCY_KEY).getAttributeAsDoubleVector();
            bool flag = false;
            for(double i : tmp) {
                if(max < i) {
                    max = i;
                    flag = true;
                }
            }
            if(flag) {
                tumorGenotype = k;
            }
        }
    }
    std::optional<std::string> phasingString = makePhasingString(tumorGenotype);
    if(!phasingString.has_value()) {
        return 0.0;
    }
    std::vector<std::pair<int, double>> phasedProbs;
    if(phasedProbabilities.find(phasingString.value()) == phasedProbabilities.end()) {
        return 0.0;
    } else {
        phasedProbs = phasedProbabilities.at(phasingString.value());
    }
    if(phasedProbs.empty()) {
        return 0.0;
    }
    max = 0.0;
    for(auto & k : phasedProbs) {
        if((k.first - k.second) <= maxIntraHaplotypeDistance) {
            if(max < k.second) {
                max = k.second;
            }
        }
    }
    return max;
}

std::optional<std::string> FilteredHaplotypeFilter::makePhasingString(std::shared_ptr<Genotype> &genotype) {
    std::string pgt;
    std::string pid;
    if(genotype->hasExtendedAttribute(VCFConstants::HAPLOTYPE_CALLER_PHASING_GT_KEY)) {
        pgt = genotype->getExtendedAttribute(VCFConstants::HAPLOTYPE_CALLER_PHASING_GT_KEY, nullptr).getAttributeAsString();
    }
    if(genotype->hasExtendedAttribute(VCFConstants::HAPLOTYPE_CALLER_PHASING_ID_KEY)) {
        pid = genotype->getExtendedAttribute(VCFConstants::HAPLOTYPE_CALLER_PHASING_ID_KEY, nullptr).getAttributeAsString();
    }
    return (pgt.size() == 0 && pid.size() == 0) ? std::optional<std::string>() : std::optional<std::string>(pgt+pid);
}

FilteredHaplotypeFilter::FilteredHaplotypeFilter(double maxIntraHaplotypeDistance) : maxIntraHaplotypeDistance(maxIntraHaplotypeDistance){

}

void FilteredHaplotypeFilter::accumulateDataForLearning(const std::shared_ptr<VariantContext> &vc,
                                                        ErrorProbabilities errorProbabilities,
                                                        Mutect2FilteringEngine *filteringEngine) {
    double artifactProbability = 0.0;
    for(auto & k : errorProbabilities.getProbabilitiesByFilter()) {
        if(k.first->errorType() != SEQUENCING) {
            if(k.first->filterName() != filterName()) {
                if(k.second > artifactProbability) {
                    artifactProbability = k.second;
                }
            }
        }
    }

    for(auto & k : *vc->getGenotypes()->getGenotypes()) {
        if(!filteringEngine->isTumor(k.get())) {
            continue;
        }
        std::optional<std::string> phasingString = makePhasingString(k);
        if(!phasingString.has_value()) {
            continue;
        }

        if(accumulatingPhasedProbabilities.find(phasingString.value()) == accumulatingPhasedProbabilities.end()) {
            accumulatingPhasedProbabilities.insert({phasingString.value(), std::vector<std::pair<int, double>>()});
        }

        accumulatingPhasedProbabilities[phasingString.value()].emplace_back(std::pair(vc->getStart(), artifactProbability));
    }
}

void FilteredHaplotypeFilter::learnParameters() {
    phasedProbabilities = accumulatingPhasedProbabilities;
}

void FilteredHaplotypeFilter::clearAccumulatedData() {
   accumulatingPhasedProbabilities.clear();
}

std::string FilteredHaplotypeFilter::filterName() {
    return VCFConstants::BAD_HAPLOTYPE_FILTER_NAME;
}

ErrorType FilteredHaplotypeFilter::errorType() {
    return ARTIFACT;
}

std::vector<std::string> FilteredHaplotypeFilter::requiredAnnotations() {
    return {};
}

int FilteredHaplotypeFilter::filterIndex() {
    return 15;
}

