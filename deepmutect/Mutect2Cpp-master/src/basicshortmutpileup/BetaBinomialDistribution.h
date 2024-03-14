//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_BETABINOMIALDISTRIBUTION_H
#define MUTECT2CPP_MASTER_BETABINOMIALDISTRIBUTION_H

#include <vector>

class BetaBinomialDistribution {
private:
    double alpha;
    double beta;
    int n;
    int rng;
    constexpr static double HALF_LOG_TWO_PI = .9189385332046727;

    static double sumDeltaMinusDeltaSum(double p, double q);
    static double deltaMinusDeltaSum(double a, double b);
    static double logGammaMinusLogGammaSum(double a, double b);
    static double logGammaSum(double a, double b);

    static std::vector<double> DELTA;

public:
    BetaBinomialDistribution(double alpha, double beta, int n, int rng);
    double logProbability(int k);
    static double logBeta(double p, double q);
};


#endif //MUTECT2CPP_MASTER_BETABINOMIALDISTRIBUTION_H
