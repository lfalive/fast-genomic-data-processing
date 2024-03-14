//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_BETABINOMIALCLUSTER_H
#define MUTECT2CPP_MASTER_BETABINOMIALCLUSTER_H

#include "AlleleFractionCluster.h"
#include "BetaDistributionShape.h"


class BetaBinomialCluster : public AlleleFractionCluster{
private:
    static const double RATE;
    static const int NUM_EPOCHS;
    BetaDistributionShape betaDistributionShape;

    static double logOddsCorrection(const BetaDistributionShape& originalBeta, const BetaDistributionShape& newBeta, int altCount, int refCount);

public:
    explicit BetaBinomialCluster(const BetaDistributionShape& betaDistributionShape);

    virtual ~BetaBinomialCluster();

    static double logLikelihood(const Datum& datum, const BetaDistributionShape& betaDistributionShape);

    double logLikelihood(const Datum& datum) override;

    double logLikelihood(int totalCount, int altCount) override;

    void learn(std::vector<Datum> data) override;

    std::string toString() override;
};


#endif //MUTECT2CPP_MASTER_BETABINOMIALCLUSTER_H
