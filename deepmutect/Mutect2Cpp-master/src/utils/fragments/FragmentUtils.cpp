//
// Created by lhh on 5/7/22.
//

#include <cassert>
#include "FragmentUtils.h"
#include "ReadUtils.h"

void FragmentUtils::adjustQualsOfOverlappingPairedFragments(
        std::pair<std::shared_ptr<SAMRecord>, std::shared_ptr<SAMRecord>> &pair, bool setConflictingToZero,
        int halfOfPcrSnvQual, int halfOfPcrIndelQual) {

    bool inOrder = pair.first->getSoftStart() < pair.second->getSoftStart();
    std::shared_ptr<SAMRecord> firstRead = inOrder ? pair.first : pair.second;
    std::shared_ptr<SAMRecord> secondRead = inOrder ? pair.second : pair.first;

    assert(firstRead != nullptr);
    assert(secondRead != nullptr);
    assert(firstRead->getName() == secondRead->getName());

    // don't adjust fragments that do not overlap
    if (firstRead->getEnd() < secondRead->getStart() || firstRead->getContig() != secondRead->getContig()) {
        return;
    }

    std::pair<int, bool> offset = ReadUtils::getReadCoordinateForReferenceCoordinate(firstRead, secondRead->getStart());
    int firstReadStop = offset.second ? offset.first+1 : offset.first;
    int numOverlappingBases = std::min(firstRead->getLength() - firstReadStop, secondRead->getLength());

    std::shared_ptr<uint8_t[]> & firstReadBases = firstRead->getBasesNoCopy();
    std::shared_ptr<uint8_t[]> & firstReadQuals = firstRead->getBaseQualitiesNoCopy();
    std::shared_ptr<uint8_t[]> & secondReadBases = secondRead->getBasesNoCopy();
    std::shared_ptr<uint8_t[]> & secondReadQuals = secondRead->getBaseQualitiesNoCopy();

    uint8_t halfOfPcrErrorQual = halfOfPcrSnvQual;
    for(int i=0; i<numOverlappingBases; i++)
    {
        int firstReadIndex = firstReadStop + i;
        uint8_t firstReadBase = firstReadBases[firstReadIndex];
        uint8_t secondReadBase = secondReadBases[i];

        if (firstReadBase == secondReadBase) {
            firstReadQuals[firstReadIndex] = std::min(firstReadQuals[firstReadIndex], halfOfPcrErrorQual);
            secondReadQuals[i] = std::min(secondReadQuals[i], halfOfPcrErrorQual);
        } else if(setConflictingToZero) {
            // If downstream processing forces read pairs to support the same haplotype, setConflictingToZero should be false
            // because the original base qualities of conflicting bases, when pegged to the same haplotype, will
            // automatically weaken the strength of one another's evidence.  Furthermore, if one base if low quality
            // and one is high it will essentially ignore the low quality base without compromising the high-quality base
            firstReadQuals[firstReadIndex] = 0;
            secondReadQuals[i] = 0;
        }
    }

    if(halfOfPcrIndelQual != MISSING_VALUE)
    {
        uint8_t maxIndelQual = halfOfPcrIndelQual;
        int firstReadInsertionQualsLength = 0;
        int firstReadDeletionQualsLength = 0;
        int secondReadInsertionQualsLength = 0;
        int secondReadDeletionQualsLength = 0;
        auto firstReadInsertionQuals = ReadUtils::getBaseInsertionQualities(firstRead, firstReadInsertionQualsLength);
        auto firstReadDeletionQuals = ReadUtils::getBaseInsertionQualities(firstRead, firstReadDeletionQualsLength);
        auto secondReadInsertionQuals = ReadUtils::getBaseInsertionQualities(secondRead, secondReadInsertionQualsLength);
        auto secondReadDeletionQuals = ReadUtils::getBaseInsertionQualities(secondRead, secondReadDeletionQualsLength);

        for(int i=0; i<numOverlappingBases; i++)
        {
            int firstReadIndex = firstReadStop + i;
            firstReadDeletionQuals[firstReadIndex] = std::min(firstReadDeletionQuals[firstReadIndex], maxIndelQual);
            firstReadInsertionQuals[firstReadIndex] = std::min(firstReadInsertionQuals[firstReadIndex], maxIndelQual);
            secondReadDeletionQuals[i] = std::min(secondReadDeletionQuals[i], maxIndelQual);
            secondReadInsertionQuals[i] = std::min(secondReadInsertionQuals[i], maxIndelQual);
        }

        ReadUtils::setDeletionBaseQualities(firstRead, firstReadDeletionQuals, firstReadDeletionQualsLength);
        ReadUtils::setInsertionBaseQualities(firstRead, firstReadInsertionQuals, firstReadInsertionQualsLength);
        ReadUtils::setDeletionBaseQualities(secondRead, secondReadDeletionQuals, secondReadDeletionQualsLength);
        ReadUtils::setInsertionBaseQualities(secondRead, secondReadInsertionQuals, secondReadInsertionQualsLength);

    }
}