//
// Created by cluster on 22-11-11.
//

#include "BinomialCluster.h"
#include "BetaBinomialDistribution.h"
#include "BetaBinomialCluster.h"
#include<algorithm>

BinomialCluster::BinomialCluster(double mean) : betaDistributionShape(getFuzzyBinomial(mean)){

}

BetaDistributionShape BinomialCluster::getFuzzyBinomial(double unboundedMean) {
    double mean = std::min(unboundedMean, 1 - STD_DEV_OVER_MEAN);
    double alphaPlusBeta = ((1 - mean) / (mean * (STD_DEV_OVER_MEAN * STD_DEV_OVER_MEAN))) - 1;
    double alpha = mean * alphaPlusBeta;
    double beta = alphaPlusBeta - alpha;
    return BetaDistributionShape(alpha, beta);
}

void BinomialCluster::learn(std::vector<Datum> data) {
    double  alt = 0.0001;
    double  total = 0.0001;
    for(const auto & k : data) {
        alt += k.getAltCount();
        total += k.getTotalCount();
    }
    betaDistributionShape = getFuzzyBinomial(alt/total);
}

double BinomialCluster::logLikelihood(int totalCount, int altCount) {
    return BetaBinomialDistribution(betaDistributionShape.getAlpha(), betaDistributionShape.getBeta(), totalCount, 0).logProbability(altCount);
}

std::string BinomialCluster::toString() {
    return "";
}

double BinomialCluster::logLikelihood(const Datum &datum) {
    return BetaBinomialCluster::logLikelihood(datum, betaDistributionShape);
}
