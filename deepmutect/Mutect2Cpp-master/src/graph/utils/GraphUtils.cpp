//
// Created by 梦想家xixi on 2021/11/19.
//

#include "GraphUtils.h"
#include "Mutect2Utils.h"

int GraphUtils::minKmerLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>>&  kmers) {
    if(kmers.empty())
        return 0;
    int ret = kmers.begin()->second;
    for(std::pair<std::shared_ptr<uint8_t[]>, int> kmer : kmers) {
        if(kmer.second < ret)
            ret = kmer.second;
    }
    return ret;
}

int GraphUtils::commonMaximumPrefixLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> & listOfBytes) {
    int minLength = minKmerLength(listOfBytes);
    for (int i = 0; i < minLength; i++) {
        std::list<std::pair<std::shared_ptr<uint8_t[]>, int>>::iterator liter = listOfBytes.begin();
        uint8_t b = listOfBytes.begin()->first.get()[i];
        liter++;
        for(; liter != listOfBytes.end(); liter++) {
            if(b != liter->first.get()[i]) {
                return i;
            }
        }
    }
    return minLength;
}

int GraphUtils::commonMaximumSuffixLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> & listOfBytes,  const int minLength) {
    for(int suffixLen = 0; suffixLen < minLength; suffixLen++) {
        uint8_t b = listOfBytes.begin()->first.get()[listOfBytes.begin()->second - suffixLen -1];
        std::list<std::pair<std::shared_ptr<uint8_t[]>, int>>::iterator liter = listOfBytes.begin();
        liter++;
        for(; liter != listOfBytes.end(); liter++) {
            if(b != liter->first.get()[liter->second - suffixLen - 1]) {
                return suffixLen;
            }
        }
    }
    return minLength;
}

std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> GraphUtils::getKmers(const std::vector<std::shared_ptr<SeqVertex>>& vertices) {
    Mutect2Utils::validateArg(!vertices.empty(), "no vertex");
    std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> ret;
    for(std::shared_ptr<SeqVertex> v : vertices) {
        ret.emplace_back(std::pair<std::shared_ptr<uint8_t[]>, int>(v->getSequence(), v->getLength()));
    }
    return ret;
}

std::list<std::pair<std::shared_ptr<uint8_t[]>, int>>
GraphUtils::getKmers(const phmap::flat_hash_set<std::shared_ptr<SeqVertex>>& vertices) {
    Mutect2Utils::validateArg(!vertices.empty(), "no vertex");
    std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> ret;
    for(std::shared_ptr<SeqVertex> v : vertices) {
        ret.emplace_back(std::pair<std::shared_ptr<uint8_t[]>, int>(v->getSequence(), v->getLength()));
    }
    return ret;
}
