//
// Created by lhh on 6/15/22.
//

#include <cassert>
#include "BaseQualityRankSumTest.h"
#include "ReadUtils.h"

std::optional<double> BaseQualityRankSumTest::getReadBaseQuality(std::shared_ptr<SAMRecord> read, int refLoc) {
    assert(read != nullptr);

    int readCoordinate = ReadUtils::getReadCoordinateForReferenceCoordinate(read, refLoc, ClippingTail::RIGHT_TAIL, true);
    return readCoordinate == ReadUtils::CLIPPING_GOAL_NOT_REACHED ? std::nullopt : std::optional<double>(read->getBaseQuality(readCoordinate));
}