//
// Created by 梦想家xixi on 2021/12/1.
//

#include "KMerCounter.h"
#include "Mutect2Utils.h"

KMerCounter::KMerCounter(const int kmerLength) : kmerLength(kmerLength){
    Mutect2Utils::validateArg(kmerLength > 0, "kmerLength must be > 0");
}

void KMerCounter::addKmer(Kmer *kmer, int kmerCount) {
    Mutect2Utils::validateArg(kmer->getLength() == kmerLength, "bad kmer length");
    Mutect2Utils::validateArg(kmerCount >= 0, "bad kmerCount");

    CountedKmer* countFromMap = countsByKMer.at(kmer);
    if(countFromMap == nullptr) {
        countFromMap = new CountedKmer(kmer);
        countsByKMer.insert(std::pair<Kmer*, CountedKmer*>(kmer, countFromMap));
    }
    countFromMap->count += kmerCount;
}
