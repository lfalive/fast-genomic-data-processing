//
// Created by 梦想家xixi on 2021/12/23.
//

#include "ReadCoordinateComparator.h"
#include "Mutect2Utils.h"
#include "ReadUtils.h"

ReadCoordinateComparator::ReadCoordinateComparator(SAMFileHeader *header) : header(header){
    Mutect2Utils::validateArg(header, "null is not allowed here");
}

int ReadCoordinateComparator::compare(std::shared_ptr<SAMRecord> first, std::shared_ptr<SAMRecord> second) {
    int result = compareCoordinates(first, second, header);
    if ( result != 0 ) {
        return result;
    }
    if (first->isReverseStrand() != second->isReverseStrand()) {
        return first->isReverseStrand()? 1: -1;
    }
    if(!first->getName().empty() && !second->getName().empty()) {
        if(first->getName() != second->getName()) {
            return first->getName() < second->getName() ? -1 : 1;
        }
    }
    result = Mutect2Utils::Int_compare(ReadUtils::getSAMFlagsForRead(first), ReadUtils::getSAMFlagsForRead(second));
    if(result != 0) {return result;}
    result = Mutect2Utils::Int_compare(first->getMappingQuality(), second->getMappingQuality());
    if(result != 0) {return result;}
    if(first->isPaired() && second->isPaired()) {
        result = Mutect2Utils::Int_compare(ReadUtils::getMateReferenceIndex(first, header), ReadUtils::getMateReferenceIndex(second, header));
        if(result != 0) {return result;}
        result = Mutect2Utils::Int_compare(first->getMateStart(), second->getMateStart());
        if(result != 0) {return result;}
    }
    result = Mutect2Utils::Int_compare(first->getFragmentLength(), second->getFragmentLength());
    return result;
}

int ReadCoordinateComparator::compareCoordinates(std::shared_ptr<SAMRecord> first, std::shared_ptr<SAMRecord> second, SAMFileHeader *header) {
    int firstRefIndex = ReadUtils::getAssignedReferenceIndex(first, header);
    int secondRefIndex = ReadUtils::getAssignedReferenceIndex(second, header);
    if ( firstRefIndex == -1 ) {
        return (secondRefIndex == -1 ? 0 : 1);
    }
    else if ( secondRefIndex == -1 ) {
        return -1;
    }

    int refIndexDifference = firstRefIndex - secondRefIndex;
    if ( refIndexDifference != 0 ) {
        return refIndexDifference;
    }

    return first->getAssignedStart() < second->getAssignedStart() ? -1 : (first->getAssignedStart() == second->getAssignedStart() ? 0 : 1);
}
