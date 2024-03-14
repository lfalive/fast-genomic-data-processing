//
// Created by 梦想家xixi on 2021/12/14.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGASSEMBLERARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_READTHREADINGASSEMBLERARGUMENTCOLLECTION_H

#include <string>
#include <vector>
#include "ReadThreadingAssembler.h"


class ReadThreadingAssemblerArgumentCollection {
public:
    static double DEFAULT_PRUNING_LOG_ODDS_THRESHOLD;
    static std::string ERROR_CORRECT_READS_LONG_NAME;
    static std::string CAPTURE_ASSEMBLY_FAILURE_BAM_LONG_NAME;
    std::vector<int> kmerSizes{10,25};
    bool dontIncreaseKmerSizesForCycles = false;
    bool allowNonUniqueKmersInRef = false;
    int numPruningSamples = 1;
    int minDanglingBranchLength = 4;
    bool recoverAllDanglingBranches = false;
    int maxNumHaplotypesInPopulation = 128;
    int minPruneFactor = 2;
    double initialErrorRateForPruning = 0.001;
    double pruningLogOddsThreshold = DEFAULT_PRUNING_LOG_ODDS_THRESHOLD;
    int maxUnprunedVariants = 100;
    bool debugAssembly;
    bool debugGraphTransformations = false;
    std::string graphOutput = "";
    bool captureAssemblyFailureBAM = false;
    bool errorCorrectReads = false;
    int kmerLengthForReadErrorCorrection = 25;
    int minObservationsForKmerToBeSolid = 20;
    virtual ReadThreadingAssembler* makeReadThreadingAssembler() = 0;
    virtual bool consensusMode() { return false; }

    int ggaExtension = 300;
    int discoverExtension = 25;
    int snpPadding = 20;
    int indelPadding = 150;
    bool dontTrimActiveRegions = false;

};


#endif //MUTECT2CPP_MASTER_READTHREADINGASSEMBLERARGUMENTCOLLECTION_H
