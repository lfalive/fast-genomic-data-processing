//
// Created by 梦想家xixi on 2021/12/14.
//

#include "ReadThreadingAssemblerArgumentCollection.h"
#include "MathUtils.h"

double ReadThreadingAssemblerArgumentCollection::DEFAULT_PRUNING_LOG_ODDS_THRESHOLD = MathUtils::log10ToLog(1.0);
std::string ReadThreadingAssemblerArgumentCollection::ERROR_CORRECT_READS_LONG_NAME("error-correct-reads");
std::string ReadThreadingAssemblerArgumentCollection::CAPTURE_ASSEMBLY_FAILURE_BAM_LONG_NAME("capture-assembly-failure-bam");
