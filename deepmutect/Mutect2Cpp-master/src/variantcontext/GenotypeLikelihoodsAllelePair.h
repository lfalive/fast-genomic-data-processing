//
// Created by 梦想家xixi on 2021/11/27.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPELIKELIHOODSALLELEPAIR_H
#define MUTECT2CPP_MASTER_GENOTYPELIKELIHOODSALLELEPAIR_H


class GenotypeLikelihoodsAllelePair {
public:
    const int alleleIndex1;
    const int alleleIndex2;
    GenotypeLikelihoodsAllelePair(int alleleIndex1, int alleleIndex2) : alleleIndex1(alleleIndex1), alleleIndex2(alleleIndex2){}
    GenotypeLikelihoodsAllelePair() : alleleIndex1(0), alleleIndex2(0) {};
};


#endif //MUTECT2CPP_MASTER_GENOTYPELIKELIHOODSALLELEPAIR_H
