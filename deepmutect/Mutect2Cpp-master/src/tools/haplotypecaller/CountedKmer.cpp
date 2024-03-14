//
// Created by 梦想家xixi on 2021/12/1.
//

#include "CountedKmer.h"

bool CountedKmer::operator<(const CountedKmer &o) const {
    return count < o.count;
}
