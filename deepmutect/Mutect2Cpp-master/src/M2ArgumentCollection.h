//
// Created by lhh on 10/23/21.
//

#ifndef MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H

#include <set>
#include "haplotypecaller/LikelihoodEngineArgumentCollection.h"

typedef struct M2ArgumentCollection{
    int callableDepth;
    int maxProbPropagationDistance;
    double activeProbThreshold;
    int assemblyRegionPadding;
    int minAssemblyRegionSize;
    int maxAssemblyRegionSize;
    //std::set<std::string> normalSamples;
    std::string normalSample;
    bool genotypeGermlineSites = false;
    LikelihoodEngineArgumentCollection likelihoodArgs;
    int pcrSnvQual = 40;
    int pcrIndelQual = 40;
    int maxMnpDistance = 1;
    double minAF = 0.0;
    bool mitochondria = false;
    double emissionLog10Odds = DEFAULT_EMISSION_LOG_10_ODDS;
    double normalLog10Odds = DEFAULT_NORMAL_LOG_10_ODDS;
    double afOfAllelesNotInGermlineResource = -1;


    constexpr static double DEFAULT_EMISSION_LOG_10_ODDS = 3.0;
    constexpr static double DEFAULT_MITO_EMISSION_LOD = 0;
    constexpr static double DEFAULT_NORMAL_LOG_10_ODDS = 2.2;
    constexpr static double DEFAULT_AF_FOR_MITO_CALLING = 4e-3;
    constexpr static double DEFAULT_AF_FOR_TUMOR_ONLY_CALLING = 5e-8;
    constexpr static double DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING = 1e-6;


    static double getInitialLogOdds() {
		// A constant should not be returned, it caused a bug sometimes, one over a billion
		// So I won't implement MTAC class for the time being
		// 4.605170185988092 == 2 * ln(10)
		return 4.605170185988092;
    }

    double getEmissionLogOdds(){
        return MathUtils::log10ToLog(mitochondria && emissionLog10Odds == DEFAULT_EMISSION_LOG_10_ODDS ? DEFAULT_MITO_EMISSION_LOD : emissionLog10Odds);
    }

    double getDefaultAlleleFrequency(){
        return afOfAllelesNotInGermlineResource >= 0 ? afOfAllelesNotInGermlineResource :
        (mitochondria ? DEFAULT_AF_FOR_MITO_CALLING:
        (normalSample.empty() ? DEFAULT_AF_FOR_TUMOR_ONLY_CALLING : DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING));
    }
}M2ArgumentCollection;

#endif //MUTECT2CPP_MASTER_M2ARGUMENTCOLLECTION_H
