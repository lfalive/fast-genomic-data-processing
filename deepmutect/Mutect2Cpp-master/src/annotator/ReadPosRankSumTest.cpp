//
// Created by lhh on 6/15/22.
//

#include <cassert>
#include "ReadPosRankSumTest.h"
#include "ReadUtils.h"
#include "AlignmentUtils.h"

std::optional<double>
ReadPosRankSumTest::getReadPosition(std::shared_ptr<SAMRecord> read, int refLoc) {
    assert(read != nullptr);
    int offset = ReadUtils::getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), refLoc, ClippingTail::RIGHT_TAIL, true);
    if ( offset == ReadUtils::CLIPPING_GOAL_NOT_REACHED ) {
        return std::nullopt;
    }

    // If the offset inside a deletion, it does not lie on a read.
    if ( AlignmentUtils::isInsideDeletion(read->getCigar(), offset) ) {
        return std::optional<double>(INVALID_ELEMENT_FROM_READ);
    }

    // hard clips at this point in the code are perfectly good bases that were clipped to make the read fit the assembly region
    auto cigar =  read->getCigar();
    CigarElement firstElement = cigar->getFirstCigarElement();
    CigarElement lastElement = cigar->getLastCigarElement();
    int leadingHardClips = firstElement.getOperator() == CigarOperator::H ? firstElement.getLength() : 0;
    int trailingHardClips = lastElement.getOperator() == CigarOperator::H ? lastElement.getLength() : 0;
    int readPos = leadingHardClips + AlignmentUtils::calcAlignmentByteArrayOffset(read->getCigar(), offset, false, 0, 0);
    int numAlignedBases = AlignmentUtils::getNumAlignedBasesCountingSoftClips( read );
    int numOriginalBases = numAlignedBases + leadingHardClips + trailingHardClips;

    //After the middle of the read, we compute the postion from the end of the read.
    if (readPos > numOriginalBases / 2) {
        readPos = numOriginalBases - (readPos + 1);
    }

    return {readPos};
}