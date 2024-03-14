//
// Created by 梦想家xixi on 2021/12/15.
//

#include "MutectReadThreadingAssemblerArgumentCollection.h"

ReadThreadingAssembler *MutectReadThreadingAssemblerArgumentCollection::makeReadThreadingAssembler() {
    ReadThreadingAssembler* assemblyEngine = new ReadThreadingAssembler(maxNumHaplotypesInPopulation, kmerSizes,
                                                                        dontIncreaseKmerSizesForCycles, allowNonUniqueKmersInRef, numPruningSamples, disableAdaptivePruning ? minPruneFactor : 0,
                                                                        !disableAdaptivePruning, initialErrorRateForPruning, pruningLogOddsThreshold, maxUnprunedVariants);

    assemblyEngine->setRecoverDanglingBranches(true);
    assemblyEngine->setRecoverAllDanglingBranches(recoverAllDanglingBranches);
    assemblyEngine->setMinDanglingBranchLength(minDanglingBranchLength);

    return assemblyEngine;
}
