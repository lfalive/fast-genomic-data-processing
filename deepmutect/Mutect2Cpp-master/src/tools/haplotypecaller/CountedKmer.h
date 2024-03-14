//
// Created by 梦想家xixi on 2021/12/1.
//

#ifndef MUTECT2CPP_MASTER_COUNTEDKMER_H
#define MUTECT2CPP_MASTER_COUNTEDKMER_H


#include "Kmer.h"

class CountedKmer {
public:
    Kmer* kmer;

    Kmer* getKmer() { return kmer;}

    int getCount() const {return count;}

    bool operator<(const CountedKmer & o) const;

    CountedKmer(Kmer* kmer) : kmer(kmer) {}

    int count = 0;
};


#endif //MUTECT2CPP_MASTER_COUNTEDKMER_H
