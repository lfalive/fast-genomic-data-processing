//
// Created by 梦想家xixi on 2022/1/11.
//

#include "ReadPileup.h"


ReadPileup::ReadPileup(int tid, int pos, const std::list<pileRead*> & reads) : tid(tid), pos(pos), reads(reads){
}

const std::list<pileRead*> & ReadPileup::getPileupElements() {
    return reads;
}

int ReadPileup::getPosition() {
    return pos;
}
