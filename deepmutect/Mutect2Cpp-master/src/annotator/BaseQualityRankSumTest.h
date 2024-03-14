//
// Created by lhh on 6/15/22.
//

#ifndef MUTECT2CPP_MASTER_BASEQUALITYRANKSUMTEST_H
#define MUTECT2CPP_MASTER_BASEQUALITYRANKSUMTEST_H

#include <optional>
#include "samtools/SAMRecord.h"

class BaseQualityRankSumTest {
public:
    static std::optional<double> getReadBaseQuality(std::shared_ptr<SAMRecord> read, int refLoc);
};


#endif //MUTECT2CPP_MASTER_BASEQUALITYRANKSUMTEST_H
