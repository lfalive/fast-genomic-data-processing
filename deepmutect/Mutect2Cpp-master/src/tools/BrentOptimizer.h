//
// Created by cluster on 22-11-9.
//

#ifndef MUTECT2CPP_MASTER_BRENTOPTIMIZER_H
#define MUTECT2CPP_MASTER_BRENTOPTIMIZER_H
#include "UnivariatePointValuePair.h"
#include <vector>
#include "EStep.h"

class BrentOptimizer {
private:
    double relativeThreshold;
    double absoluteThreshold;
    double min;
    double max;
    double start;
    double (*fun)(double, double, const std::vector<EStep> &);
    int maxEval;
    int currentEval;
    double artifactBetaMean;
    std::vector<EStep> potentialArtifacts;
    static const double GOLDEN_SECTION;

    double computeObjectiveValue(double x);
    UnivariatePointValuePair best(UnivariatePointValuePair* a, UnivariatePointValuePair* b, bool isMinim);

public:
    BrentOptimizer(double rel, double abs, double min, double max, double guess, double (*fun)(double, double, const std::vector<EStep> &), int maxEval, double artifactBetaMean, const std::vector<EStep> & potentialArtifacts);
    UnivariatePointValuePair doOptimize();
    bool AreDoubleSame(double dFirstVal, double dSecondVal);
};


#endif //MUTECT2CPP_MASTER_BRENTOPTIMIZER_H
