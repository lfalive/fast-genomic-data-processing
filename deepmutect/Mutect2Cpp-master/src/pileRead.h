//
// Created by 梦想家xixi on 2022/2/24.
//

#ifndef MUTECT2CPP_MASTER_PILEREAD_H
#define MUTECT2CPP_MASTER_PILEREAD_H

#include "samtools/SAMRecord.h"

typedef struct {
    std::shared_ptr<SAMRecord> read;
    int activateStart;  // When traversing, the position greater than (or equal to) this variable will be considered, otherwise not considered
    int activateStop;
}pileRead;

#endif //MUTECT2CPP_MASTER_PILEREAD_H
