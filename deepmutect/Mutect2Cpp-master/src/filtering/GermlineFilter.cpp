//
// Created by cluster on 22-11-14.
//

#include "GermlineFilter.h"
#include "MathUtils.h"
#include "BinomialDistribution.h"
#include "NaturalLogUtils.h"

double GermlineFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                 Mutect2FilteringEngine *filteringEngine,
                                                 std::shared_ptr<ReferenceContext>) {
    std::vector<double> somaticLogOdds = Mutect2FilteringEngine::getTumorLogOdds(vc);
    int maxLodIndex = MathUtils::maxElementIndex(somaticLogOdds);
    std::vector<double> tmp;
    if(vc->hasAttribute(VCFConstants::NORMAL_LOG_10_ODDS_KEY)) {
       tmp = vc->getAttributes().at(VCFConstants::NORMAL_LOG_10_ODDS_KEY).getAttributeAsDoubleVector();
    }
    MathUtils::applyToArrayInPlace(tmp, MathUtils::log10ToLog);
    std::optional<std::vector<double>> normalLogOdds = vc->hasAttribute(VCFConstants::NORMAL_LOG_10_ODDS_KEY) ? tmp : std::optional<std::vector<double>>();
    std::vector<double> negativeLog10AlleleFrequencies;
    if(vc->hasAttribute(VCFConstants::POPULATION_AF_KEY)) {
        negativeLog10AlleleFrequencies = vc->getAttributes().at(VCFConstants::POPULATION_AF_KEY).getAttributeAsDoubleVector();
    } else {
        negativeLog10AlleleFrequencies = {-1};
    }
    double populationAF = std::pow(10, -negativeLog10AlleleFrequencies[maxLodIndex]);

    if (populationAF < EPSILON) {
        return 0;
    } else if (populationAF > 1 - EPSILON) {
        return 1;
    }

    std::vector<int> alleleCounts = filteringEngine->sumADsOverSamples(vc, true, false);
    int totalCount = 0;
    for(int i : alleleCounts) {
        totalCount += i;
    }
    if (totalCount == 0) {  // this could come up in GGA mode
        return 0;
    }
    int refCount = alleleCounts[0];
    int altCount = alleleCounts[maxLodIndex + 1];
    double altAlleleFraction = filteringEngine->weightedAverageOfTumorAFs(vc)[maxLodIndex];

    double maf = computeMinorAlleleFraction(vc, filteringEngine, alleleCounts);
    tmp = {BinomialDistribution(totalCount, maf).logProbability(altCount), BinomialDistribution(totalCount, 1- maf).logProbability(altCount)};
    double logGermlineLikelihood = NaturalLogUtils::LOG_ONE_HALF + NaturalLogUtils::logSumExp(tmp);
    double logSomaticLikelihood = filteringEngine->getSomaticClusteringModel()->logLikelihoodGivenSomatic(totalCount, altCount);
    double logOddsOfGermlineHetVsSomatic = logGermlineLikelihood - logSomaticLikelihood;
    double logOddsOfGermlineHomAltVsSomatic = altAlleleFraction < MIN_ALLELE_FRACTION_FOR_GERMLINE_HOM_ALT ? (-1.0 / 0.0) : 0;
    double normalLod = normalLogOdds.has_value() ? normalLogOdds.value()[maxLodIndex] : 0;
    return germlineProbability(-normalLod, logOddsOfGermlineHetVsSomatic, logOddsOfGermlineHomAltVsSomatic,
                               populationAF, filteringEngine->getLogSomaticPrior(vc, maxLodIndex));
}

double GermlineFilter::computeMinorAlleleFraction(const shared_ptr<VariantContext> &vc,
                                                  Mutect2FilteringEngine *filteringEngine,
                                                  const vector<int> &alleleCounts) {
    return 0.5;
}

double GermlineFilter::germlineProbability(double normalLogOdds, double logOddsOfGermlineHetVsSomatic,
                                           double logOddsOfGermlineHomAltVsSomatic, double populationAF,
                                           double logPriorSomatic) {
    double logPriorNotSomatic = NaturalLogUtils::log1mexp(logPriorSomatic);
    double logPriorGermlineHet = std::log(2*populationAF*(1-populationAF));
    double logPriorGermlineHomAlt = std::log( populationAF * populationAF);
    double logPriorNotGermline = std::log((1 - populationAF) * (1 - populationAF));

    // the following are unnormalized probabilities
    double logProbGermlineHet = logPriorGermlineHet + logOddsOfGermlineHetVsSomatic + normalLogOdds + logPriorNotSomatic;
    double logProbGermlineHomAlt = logPriorGermlineHomAlt + logOddsOfGermlineHomAltVsSomatic + normalLogOdds + logPriorNotSomatic;
    std::vector<double> tmp = {logProbGermlineHet, logProbGermlineHomAlt};
    double logProbGermline = NaturalLogUtils::logSumExp(tmp);
    double logProbSomatic = logPriorNotGermline + logPriorSomatic;

    return NaturalLogUtils::normalizeLog({logProbGermline, logProbSomatic}, false, true)[0];
}

std::string GermlineFilter::filterName() {
    return VCFConstants::GERMLINE_QUAL_KEY;
}

std::vector<std::string> GermlineFilter::requiredAnnotations() {
    return {VCFConstants::TUMOR_LOG_10_ODDS_KEY, VCFConstants::POPULATION_AF_KEY};
}

ErrorType GermlineFilter::errorType() {
    return NON_SOMATIC;
}

int GermlineFilter::filterIndex() {
    return 4;
}
