//
// Created by cluster on 22-11-10.
//

#include "ThresholdCalculator.h"
#include <algorithm>

void ThresholdCalculator::addArtifactProbability(double artifactProbability) {
    artifactProbabilities.emplace_back(artifactProbability);
}

ThresholdCalculator::ThresholdCalculator(double initialThreshold, double maxFalseDiscoveryRate, double fScoreBeta) : threshold(initialThreshold), maxFalseDiscoveryRate(maxFalseDiscoveryRate), fScoreBeta(fScoreBeta){

}

void ThresholdCalculator::relearnThresholdAndClearAcumulatedProbabilities() {
    threshold = calculateThresholdBasedOnOptimalFScore(artifactProbabilities, fScoreBeta);;
}

double ThresholdCalculator::calculateThresholdBasedOnOptimalFScore(std::vector<double> posteriors, double beta) {
    std::sort(posteriors.begin(), posteriors.end());
    double expectedTruePositives = 0;
    for(int i = 0; i < posteriors.size(); i++) {
        expectedTruePositives += (1-posteriors[i]);
    }
    double truePositives = 0;
    double falsePositives = 0;
    double falseNegatives = expectedTruePositives;
    int optimalIndexInclusive = -1; // include all indices up to and including this. -1 mean filter all.
    double optimalFScore = 0;
    int N = posteriors.size();

    for (int n = 0; n < N; n++){
        truePositives += (1 - posteriors[n]);
        falsePositives += posteriors[n];
        falseNegatives -= 1 - posteriors[n];
        double F = (1+beta*beta)*truePositives /
                         ((1+beta*beta)*truePositives + beta*beta*falseNegatives + falsePositives);
        if (F >= optimalFScore) {
            optimalIndexInclusive = n;
            optimalFScore = F;
        }
    }

    return optimalIndexInclusive == -1 ? 0 : (optimalIndexInclusive == N - 1 ? 1 : posteriors[optimalIndexInclusive]);
}

double ThresholdCalculator::getThredshold() {
    return threshold;
}
