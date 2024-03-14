//
// Created by cluster on 22-11-8.
//

#include "StrandArtifactFilter.h"
#include "BetaBinomialDistribution.h"
#include "CombinatoricsUtils.h"
#include "MathUtils.h"
#include "tools/BrentOptimizer.h"

double StrandArtifactFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                       Mutect2FilteringEngine* filteringEngine,
                                                       std::shared_ptr<ReferenceContext>) {
    EStep probabilities = calculateArtifactProbabilities(vc, filteringEngine);
    return probabilities.getForwardArtifactResponsibility() + probabilities.getReverseArtifactResponsibility();
}

EStep StrandArtifactFilter::calculateArtifactProbabilities(const std::shared_ptr<VariantContext> &vc,
                                                           Mutect2FilteringEngine* filteringEngine) {
    std::vector<int> counts = filteringEngine->sumStrandCountsOverSamples(vc, true, false);

    int indelSize = std::abs(vc->getReference()->getLength() - vc->getAlternateAllele(0)->getLength());
    if (counts[2] + counts[3] == 0 || indelSize > LONGEST_STRAND_ARTIFACT_INDEL_SIZE) {
        return EStep(0, 0, counts[0] + counts[2], counts[1] + counts[3], counts[2], counts[3]);
    }


    return strandArtifactProbability(strandArtifactPrior, counts[0] + counts[2], counts[1] + counts[3], counts[2], counts[3], indelSize);
}

StrandArtifactFilter::StrandArtifactFilter() {
    strandArtifactPrior = INITIAL_STRAND_ARTIFACT_PRIOR;
    INITIAL_ALPHA_STRAND = 1.0;
    INITIAL_BETA_STRAND = 20.0;

    alphaStrand = INITIAL_ALPHA_STRAND;
    betaStrand = INITIAL_BETA_STRAND;
}

EStep StrandArtifactFilter::strandArtifactProbability(double strandArtifactPrior, int forwardCount, int reverseCount,
                                                      int forwardAltCount, int reverseAltCount, int indelSize) {
    double forwardLogLikelihood = artifactStrandLogLikelihood(forwardCount, forwardAltCount)
                                        + nonArtifactStrandLogLikelihood(reverseCount, reverseAltCount, indelSize);
    double reverseLogLikelihood = artifactStrandLogLikelihood(reverseCount, reverseAltCount)
                                        + nonArtifactStrandLogLikelihood(forwardCount, forwardAltCount, indelSize);
    double noneLogLikelihood = CombinatoricsUtils::binomialCoefficientLog(forwardCount, forwardAltCount)
                                     + CombinatoricsUtils::binomialCoefficientLog(reverseCount, reverseAltCount)
                                     - CombinatoricsUtils::binomialCoefficientLog(forwardCount + reverseCount, forwardAltCount + reverseAltCount)
    + BetaBinomialDistribution(1, 1, forwardCount + reverseCount, 0).logProbability(forwardAltCount + reverseAltCount);

    double forwardLogPrior = std::log(strandArtifactPrior/2);
    double reverseLogPrior = std::log(strandArtifactPrior/2);
    double noneLogPrior = std::log(1 - strandArtifactPrior);
    std::vector<double> tmp = {(forwardLogLikelihood + forwardLogPrior) * MathUtils::LOG10_OF_E,
                               (reverseLogLikelihood + reverseLogPrior)* MathUtils::LOG10_OF_E, (noneLogLikelihood + noneLogPrior) * MathUtils::LOG10_OF_E};

    std::vector<double> forwardReverseNoneProbs = MathUtils::normalizeLog10(tmp, false, true);

    return {forwardReverseNoneProbs[0], forwardReverseNoneProbs[1], forwardCount, reverseCount, forwardAltCount, reverseAltCount};
}

double StrandArtifactFilter::artifactStrandLogLikelihood(int strandCount, int strandAltCount) {
    return artifactStrandLogLikelihood(strandCount, strandAltCount, alphaStrand, betaStrand);
}

double
StrandArtifactFilter::artifactStrandLogLikelihood(int strandCount, int strandAltCount, double alpha, double beta) {
    return BetaBinomialDistribution(alpha, beta, strandCount, 0).logProbability(strandAltCount);
}

double StrandArtifactFilter::nonArtifactStrandLogLikelihood(int strandCount, int strandAltCount, int indelSize) {
    double betaSeq = indelSize == 0 ? BETA_SEQ_SNV : (indelSize < LONG_INDEL_SIZE ? BETA_SEQ_SHORT_INDEL : BETA_SEQ_LONG_INDEL);
    return BetaBinomialDistribution(ALPHA_SEQ, betaSeq, strandCount, 0).logProbability(strandAltCount);
}

std::vector<std::string> StrandArtifactFilter::requiredAnnotations() {
    return {};
}

void StrandArtifactFilter::learnParameters() {
    std::vector<EStep> potentialArtifacts;
    for(auto &es : eSteps) {
        if(es.getArtifactProbability() > 0.1) {
            potentialArtifacts.emplace_back(es);
        }
    }
    double totalArtifacts = 0.0;
    double totalNonArtifacts = 0.0;
    double artifactAltCount = 0.0;
    double artifactDepth = 0.0;
    for(auto &es : potentialArtifacts) {
        totalArtifacts += es.getArtifactProbability();
        totalNonArtifacts += (1 - es.getArtifactProbability());
        artifactAltCount += (es.getForwardArtifactResponsibility() * es.getForwardAltCount() + es.getReverseArtifactResponsibility() * es.getReverseAltCount());
        artifactDepth += (es.getForwardArtifactResponsibility() * es.getForwardCount() + es.getReverseArtifactResponsibility() * es.getReverseCount());
    }
    double artifactBetaMean = (artifactAltCount + INITIAL_ALPHA_STRAND) / (artifactDepth + INITIAL_ALPHA_STRAND + INITIAL_BETA_STRAND);
    BrentOptimizer optimizer(0.01, 0.01, 0.01, 100, INITIAL_ALPHA_STRAND, [](double x, double artifactBetaMean, const std::vector<EStep> & potentialArtifacts) {
        double beta = (1 / artifactBetaMean - 1) * x;
        double res = 0;
        for(auto & es : potentialArtifacts) {
            res += (es.getForwardArtifactResponsibility() * artifactStrandLogLikelihood(es.getForwardCount(), es.getForwardAltCount(), x, beta)
                    + es.getReverseArtifactResponsibility() * artifactStrandLogLikelihood(es.getReverseCount(), es.getReverseAltCount(), x, beta));
        }
        return res;
    }, 100, artifactBetaMean, potentialArtifacts);

    alphaStrand = optimizer.doOptimize().getPoint();
    betaStrand = (1/artifactBetaMean - 1)*alphaStrand;
    // free up memory
    eSteps.clear();
}

void StrandArtifactFilter::clearAccumulatedData() {
    eSteps.clear();
}

void StrandArtifactFilter::accumulateDataForLearning(const shared_ptr<VariantContext> &vc,
                                                     ErrorProbabilities errorProbabilities,
                                                     Mutect2FilteringEngine* filteringEngine) {
    EStep eStep = calculateArtifactProbabilities(vc, filteringEngine);
    eSteps.emplace_back(eStep);
}

ErrorType StrandArtifactFilter::errorType() {
    return ARTIFACT;
}

std::string StrandArtifactFilter::filterName() {
    return VCFConstants::STRAND_ARTIFACT_FILTER_NAME;
}

int StrandArtifactFilter::filterIndex() {
    return 3;
}

