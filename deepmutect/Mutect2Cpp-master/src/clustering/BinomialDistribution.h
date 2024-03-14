//
// Created by cluster on 22-11-14.
//

#ifndef MUTECT2CPP_MASTER_BINOMIALDISTRIBUTION_H
#define MUTECT2CPP_MASTER_BINOMIALDISTRIBUTION_H

#include <vector>

class BinomialDistribution {
private:
    int numberOfTrials;
    double probabilityOfSuccess;
    constexpr static double DEFAULT_EPSILON = 1e-14;
    static double getB(int n, double x, double a, double b);
    static double getA(int n, double x) {
        return 1.0;
    }
    double logBinomialProbability(int x, int n, double p, double q);
    static double getDeviancePart(double x, double mu);
    static double getStirlingError(double z);
    const static std::vector<double> EXACT_STIRLING_ERRORS;

public:
    BinomialDistribution(int trails, double p);
    double cumulativeProbability(int x);
    static double regularizedBeta(double x, double a, double b);
    static double regularizedBeta(double x,
                                  double a, double b,
                                  double epsilon, int maxIterations);
    static double evaluate(double x, double epsilon, int maxIterations, double _a, double _b);
    double logProbability(int x);

    double probability(int i);
};


#endif //MUTECT2CPP_MASTER_BINOMIALDISTRIBUTION_H
