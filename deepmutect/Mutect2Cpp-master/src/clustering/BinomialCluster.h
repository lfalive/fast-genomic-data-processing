//
// Created by cluster on 22-11-11.
//

#ifndef MUTECT2CPP_MASTER_BINOMIALCLUSTER_H
#define MUTECT2CPP_MASTER_BINOMIALCLUSTER_H

#include "AlleleFractionCluster.h"
#include "BetaDistributionShape.h"


class BinomialCluster : public AlleleFractionCluster{
private:
    constexpr static double STD_DEV_OVER_MEAN = 0.01;
    BetaDistributionShape betaDistributionShape;
    BetaDistributionShape getFuzzyBinomial(double unboundedMean);

public:
    BinomialCluster(double mean);

    virtual void learn(std::vector<Datum> data) override;

    virtual double logLikelihood(int totalCount, int altCount) override;

    virtual std::string toString() override;

    virtual double logLikelihood(const Datum& datum) override;
};


#endif //MUTECT2CPP_MASTER_BINOMIALCLUSTER_H
