//
// Created by 梦想家xixi on 2021/12/21.
//

#ifndef MUTECT2CPP_MASTER_CLIPPINGOP_H
#define MUTECT2CPP_MASTER_CLIPPINGOP_H

#include "samtools/SAMRecord.h"
#include "ReadClipper.h"
#include "ClippingRepresentation.h"


class ClippingOp {
public:
    int start;
    int stop;
    ClippingOp(int start, int stop);
    std::shared_ptr<SAMRecord> apply(ClippingRepresentation algorithm, std::shared_ptr<SAMRecord> originalRead, bool runAsserts);

private:
    std::shared_ptr<SAMRecord> applyHardClipBases(std::shared_ptr<SAMRecord> read, int start, int stop);
    std::shared_ptr<SAMRecord> applyRevertSoftClippedBases(const std::shared_ptr<SAMRecord>& read);
    class CigarShift {
    public:
        CigarShift(std::shared_ptr<Cigar> cigar, int shiftFromStart, int shiftFromEnd);

        int shiftFromStart;
        int shiftFromEnd;
        std::shared_ptr<Cigar>  cigar;
    };
    std::shared_ptr<CigarShift> hardClipCigar(std::shared_ptr<Cigar>  cigar, int start, int stop);
    int calculateHardClippingAlignmentShift(CigarElement & cigarElement, int clippedLength);
    std::shared_ptr<CigarShift> cleanHardClippedCigar(std::shared_ptr<Cigar> cigar);
    static int calculateAlignmentStartShift(std::shared_ptr<Cigar> oldCigar, std::shared_ptr<Cigar> newCigar);
    static int calcHardSoftOffset(std::shared_ptr<Cigar> cigar);
};



#endif //MUTECT2CPP_MASTER_CLIPPINGOP_H
