//
// Created by 梦想家xixi on 2021/11/12.
//

#ifndef MUTECT2CPP_MASTER_SMITHWATERMANALIGNMENT_H
#define MUTECT2CPP_MASTER_SMITHWATERMANALIGNMENT_H

#include <memory>
#include "cigar/Cigar.h"

class SmithWatermanAlignment {
public:
    virtual std::shared_ptr<Cigar> & getCigar() = 0;
    virtual int getAlignmentOffset() = 0;
    virtual ~SmithWatermanAlignment() {};
};


#endif //MUTECT2CPP_MASTER_SMITHWATERMANALIGNMENT_H
