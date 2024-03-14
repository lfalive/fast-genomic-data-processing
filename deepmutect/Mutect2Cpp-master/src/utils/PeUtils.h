//
// Created by 梦想家xixi on 2022/1/4.
//

#ifndef MUTECT2CPP_MASTER_PEUTILS_H
#define MUTECT2CPP_MASTER_PEUTILS_H

#include "htslib/sam.h"
#include "cigar/Cigar.h"
#include "samtools/SAMRecord.h"

class PeUtils {
public:
    static const char  DELETION_QUAL = 16;
    bool isBeforeSoftClip();
    bool isImmediatelyBefore(CigarOperator cigarOperator);
    bool isImmediatelyAfter(CigarOperator op);
    bool isAfterSoftClip();
    PeUtils(SAMRecord* pe, int pos);
    bool isDeletion();
    bool isBeforeDeletionStart();
    bool isBeforeInsertion();
    bool isOnGenomeCigar(CigarOperator cigarOperator);
    bool atEndOfCurrentCigar();
    int getLengthOfImmediatelyFollowingIndel();
    CigarElement * getCurrentCigarElement();
    CigarElement * getNearestOnGenomeCigarElement(int direction);
    CigarElement *getNextOnGenomeCigarElement();
    std::vector<CigarElement*> * getBetweenNextPosition();
    uint8_t getQual();
    uint8_t getBase();
    uint8_t getBaseQuality(int pos);

    int getOffset() const {return offset;}

private:
    int Cigar_offset;
    int pos;
    int currentStart;
    int offset;
    std::vector<CigarElement> * nCigarElements;
    CigarElement * currentCigarElement;
    SAMRecord *pe;
    CigarElement* getNextIndelCigarElement();
    std::vector<CigarElement*>* getBetween(int direction);
};


#endif //MUTECT2CPP_MASTER_PEUTILS_H
