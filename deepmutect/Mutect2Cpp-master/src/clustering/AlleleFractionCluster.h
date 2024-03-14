//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_ALLELEFRACTIONCLUSTER_H
#define MUTECT2CPP_MASTER_ALLELEFRACTIONCLUSTER_H


#include "Datum.h"
#include <vector>
#include <string>

class AlleleFractionCluster {
public:

    virtual double logLikelihood(int totalCount, int altCount) = 0;

    virtual void learn(std::vector<Datum> data) = 0;

    virtual std::string toString() = 0;

    virtual double logLikelihood(const Datum& datum) = 0;

    virtual ~AlleleFractionCluster() = default;
};


#endif //MUTECT2CPP_MASTER_ALLELEFRACTIONCLUSTER_H
