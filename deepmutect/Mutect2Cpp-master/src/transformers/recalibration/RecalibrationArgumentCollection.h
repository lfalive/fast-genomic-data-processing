//
// Created by lhh on 4/6/22.
//

/**
 * A collection of the arguments that are used for BQSR. Used to be common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 */

#ifndef MUTECT2CPP_MASTER_RECALIBRATIONARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_RECALIBRATIONARGUMENTCOLLECTION_H

#include <string>

class RecalibrationArgumentCollection {
public:
    int MISMATCHES_CONTEXT_SIZE = 2;
    int INDELS_CONTEXT_SIZE = 3;
    int MAXIMUM_CYCLE_VALUE = 500;
    char MISMATCHES_DEFAULT_QUALITY = -1;
    char INSERTIONS_DEFAULT_QUALITY = 45;
    char DELETIONS_DEFAULT_QUALITY = 45;
    char LOW_QUAL_TAIL = 2;
    int QUANTIZING_LEVELS = 16;
    std::string BINARY_TAG_NAME;

};


#endif //MUTECT2CPP_MASTER_RECALIBRATIONARGUMENTCOLLECTION_H
