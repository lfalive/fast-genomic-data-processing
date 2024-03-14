//
// Created by cluster on 22-11-15.
//

#include "PolymeraseSlippageFilter.h"
#include "BinomialDistribution.h"

PolymeraseSlippageFilter::PolymeraseSlippageFilter(int minSlippageLength, double slippageRate) : minSlippageLength(
        minSlippageLength), slippageRate(slippageRate) {

}

ErrorType PolymeraseSlippageFilter::errorType() {
    return ARTIFACT;
}

double PolymeraseSlippageFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                           Mutect2FilteringEngine *filteringEngine,
                                                           std::shared_ptr<ReferenceContext>) {
    std::vector<int> rpa = vc->getAttributeAsIntVector(VCFConstants::REPEATS_PER_ALLELE_KEY, {0});
    if (rpa.size() < 2) {
        return 0;
    }
    std::string ru;
    if (vc->hasAttribute(VCFConstants::REPEAT_UNIT_KEY)) {
        ru = vc->getAttributes().at(VCFConstants::REPEAT_UNIT_KEY).getAttributeAsString();
    }

    int referenceSTRBaseCount = ru.size() * rpa[0];
    int numPCRSlips = rpa[0] - rpa[1];
    if (referenceSTRBaseCount >= minSlippageLength && std::abs(numPCRSlips) == 1) {
        // calculate the p-value that out of n reads we would have at least k slippage reads
        // if this p-value is small we keep the variant (reject the PCR slippage hypothesis)
        std::vector<int> ADs = filteringEngine->sumADsOverSamples(vc, true, false);
        if ((int) ADs.size() < 2) {
            return 0;
        }
        int depth = 0;
        for (int AD: ADs) {
            depth += AD;
        }
        int altCount = depth - ADs[0];
        double logSomaticLikelihood = filteringEngine->getSomaticClusteringModel()->logLikelihoodGivenSomatic(depth,
                                                                                                              altCount);
        double likelihoodGivenSlippageArtifact;
        try {
            likelihoodGivenSlippageArtifact = BinomialDistribution::regularizedBeta(slippageRate, ADs[1] + 1,
                                                                                    ADs[0] + 1);
        } catch (...) {
            //if the special function can't be computed, use a binomial with fixed probability
            likelihoodGivenSlippageArtifact = BinomialDistribution(depth, slippageRate).probability(ADs[1]);
        }

        double logOdds = logSomaticLikelihood - std::log(likelihoodGivenSlippageArtifact);
        return filteringEngine->posteriorProbabilityOfError(vc, logOdds, 0);
    } else {
        return 0;
    }
}

std::string PolymeraseSlippageFilter::filterName() {
    return VCFConstants::POLYMERASE_SLIPPAGE;
}

std::vector<std::string> PolymeraseSlippageFilter::requiredAnnotations() {
    return {VCFConstants::REPEATS_PER_ALLELE_KEY, VCFConstants::REPEAT_UNIT_KEY};
}

int PolymeraseSlippageFilter::filterIndex() {
    return 8;
}
