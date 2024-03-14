//
// Created by cluster on 22-11-8.
//

#include "BetaBinomialCluster.h"
#include "Mutect2/SomaticLikelihoodsEngine.h"
#include "basicshortmutpileup/BetaBinomialDistribution.h"
#include "DigammaCache.h"
#include <sstream>

const double BetaBinomialCluster::RATE = 0.01;
const int BetaBinomialCluster::NUM_EPOCHS = 10;

BetaBinomialCluster::BetaBinomialCluster(const BetaDistributionShape &betaDistributionShape) : betaDistributionShape(betaDistributionShape){

}

double BetaBinomialCluster::logLikelihood(const Datum &datum, const BetaDistributionShape &betaDistributionShape) {
    int altCount = datum.getAltCount();
    int refCount = datum.getTotalCount() - altCount;
    return datum.getTumorLogOdds() + logOddsCorrection(BetaDistributionShape::FLAT_BETA, betaDistributionShape, altCount, refCount);
}

double
BetaBinomialCluster::logOddsCorrection(const BetaDistributionShape &originalBeta, const BetaDistributionShape &newBeta,
                                       int altCount, int refCount) {
    double res = 0;
    std::vector<double> inputs = {newBeta.getAlpha(), newBeta.getBeta()};
    res += SomaticLikelihoodsEngine::logDirichletNormalization(inputs);
    inputs = {newBeta.getAlpha() + altCount, newBeta.getBeta() + refCount};
    res -= SomaticLikelihoodsEngine::logDirichletNormalization(inputs);
    inputs = {originalBeta.getAlpha(), originalBeta.getBeta()};
    res -= SomaticLikelihoodsEngine::logDirichletNormalization(inputs);
    inputs = {originalBeta.getAlpha() + altCount, originalBeta.getBeta() + refCount};
    res += SomaticLikelihoodsEngine::logDirichletNormalization(inputs);
    return res;
}

double BetaBinomialCluster::logLikelihood(int totalCount, int altCount) {
    return BetaBinomialDistribution(betaDistributionShape.getAlpha(), betaDistributionShape.getBeta(), totalCount, 0).logProbability(altCount);
}

double BetaBinomialCluster::logLikelihood(const Datum& datum) {
    return logLikelihood(datum, betaDistributionShape);
}

void BetaBinomialCluster::learn(std::vector<Datum> data) {
    double alpha = betaDistributionShape.getAlpha();
    double beta = betaDistributionShape.getBeta();

    for (int epoch = 0; epoch < NUM_EPOCHS; epoch++) {
        for (Datum datum : data) {
            int alt = datum.getAltCount();
            int ref = datum.getTotalCount() - alt;

            double digammaOfTotalPlusAlphaPlusBeta = DigammaCache::digamma(datum.getTotalCount() + alpha + beta);
            double digammaOfAlphaPlusBeta = DigammaCache::digamma(alpha + beta);
            double alphaGradient = DigammaCache::digamma(alpha + alt) - digammaOfTotalPlusAlphaPlusBeta - DigammaCache::digamma(alpha) + digammaOfAlphaPlusBeta;
            double betaGradient = DigammaCache::digamma(beta + ref) - digammaOfTotalPlusAlphaPlusBeta - DigammaCache::digamma(beta) + digammaOfAlphaPlusBeta;

            alpha = std::max(alpha + RATE * alphaGradient, 0.5);
            beta = std::max(beta + RATE * betaGradient, 0.5);
        }
    }

    betaDistributionShape = BetaDistributionShape(alpha, beta);
}

std::string BetaBinomialCluster::toString() {
    std::ostringstream buffer;
    buffer << "alpha =" << betaDistributionShape.getAlpha() << "beta = " << betaDistributionShape.getBeta();
    return buffer.str();
}

BetaBinomialCluster::~BetaBinomialCluster() = default;
