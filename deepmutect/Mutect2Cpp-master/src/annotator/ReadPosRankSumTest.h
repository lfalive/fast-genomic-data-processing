//
// Created by lhh on 6/15/22.
//

#ifndef MUTECT2CPP_MASTER_READPOSRANKSUMTEST_H
#define MUTECT2CPP_MASTER_READPOSRANKSUMTEST_H

#include <optional>
#include <VariantContext.h>
#include <limits>
#include "samtools/SAMRecord.h"


class ReadPosRankSumTest {
protected:
    constexpr static double INVALID_ELEMENT_FROM_READ = -std::numeric_limits<double>::infinity();

public:
    static std::optional<double>  getReadPosition(std::shared_ptr<SAMRecord> read, int refLoc);
};


#endif //MUTECT2CPP_MASTER_READPOSRANKSUMTEST_H
