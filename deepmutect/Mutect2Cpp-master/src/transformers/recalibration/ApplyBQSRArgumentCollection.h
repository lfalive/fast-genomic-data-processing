//
// Created by lhh on 4/8/22.
//

#ifndef MUTECT2CPP_MASTER_APPLYBQSRARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_APPLYBQSRARGUMENTCOLLECTION_H

#include "QualityUtils.h"

/**
 * The collection of all arguments needed for ApplyBQSR.
 */
class ApplyBQSRArgumentCollection {
public:
    int PRESERVE_QSCORES_LESS_THAN = QualityUtils::MIN_USABLE_Q_SCORE;
    bool useOriginalBaseQualities;
    int quantizationLevels = 0;
    bool roundDown;
    bool emitOriginalQuals;
    double globalQScorePrior;

    ApplyBQSRArgumentCollection();
};


#endif //MUTECT2CPP_MASTER_APPLYBQSRARGUMENTCOLLECTION_H
