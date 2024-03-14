//
// Created by 梦想家xixi on 2021/12/1.
//

#ifndef MUTECT2CPP_MASTER_KMERCOUNTER_H
#define MUTECT2CPP_MASTER_KMERCOUNTER_H

#include <map>
#include "Kmer.h"
#include "CountedKmer.h"

class KMerCounter {
private:
    std::map<Kmer*, CountedKmer*> countsByKMer;
    int kmerLength;

public:
    KMerCounter(int kmerLength);
    void addKmer(Kmer* kmer, int kmerCount);
};


#endif //MUTECT2CPP_MASTER_KMERCOUNTER_H
