//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_HAPLOTYPEDATAHOLDER_H
#define MUTECT2CPP_MASTER_HAPLOTYPEDATAHOLDER_H

#include <cstdint>
#include <memory>
#include <utility>
#include <cassert>

struct HaplotypeDataHolder {
    uint8_t* haplotypeBases;
    unsigned length;
    HaplotypeDataHolder(uint8_t*  haplotypeBases, unsigned length) : haplotypeBases(haplotypeBases), length(length) {
        assert(haplotypeBases != nullptr);
    };
};


#endif //MUTECT2CPP_MASTER_HAPLOTYPEDATAHOLDER_H
