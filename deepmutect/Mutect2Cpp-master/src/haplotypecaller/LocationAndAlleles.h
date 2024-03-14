//
// Created by lhh on 5/23/22.
//

#ifndef MUTECT2CPP_MASTER_LOCATIONANDALLELES_H
#define MUTECT2CPP_MASTER_LOCATIONANDALLELES_H

#include <iostream>
#include <cassert>
#include <vector>
#include "Allele.h"

/**
 * This class exists to allow VariantContext objects to be compared based only on their location and set of alleles,
 * providing a more liberal equals method so that VariantContext objects can be placed into a Set
 * which retains only VCs that have non-redundant location and Allele lists.
 */
class LocationAndAlleles {
private:
    int loc;
    std::vector<std::shared_ptr<Allele>>& alleles;

public:
    LocationAndAlleles(int loc, std::vector<std::shared_ptr<Allele>> & alleles);

    int getLoc();

    std::vector<std::shared_ptr<Allele>> & getAlleles();
};

struct hash_LocationAndAlleles{
    size_t operator()(const std::shared_ptr<LocationAndAlleles> &locationAndAlleles) const {
        auto& alleles = locationAndAlleles->getAlleles();
        size_t allelesHashCode = 0;
        if(!alleles.empty())
        {
            for(auto& temp : alleles)
            {
                allelesHashCode += temp->hashcode();
            }
        }
        return 31 * locationAndAlleles->getLoc() + allelesHashCode;
    }
};

struct equal_LocationAndAlleles {
    bool operator()(const std::shared_ptr<LocationAndAlleles> &left, const std::shared_ptr<LocationAndAlleles> &right) const {
        assert(left != nullptr && right != nullptr);
        if(left->getLoc() != right->getLoc())
            return false;
        if(left->getAlleles().size() != right->getAlleles().size())
            return false;

        int size = left->getAlleles().size();
        auto& leftAlleles = left->getAlleles();
        auto& rightAlleles = right->getAlleles();
        for(int i=0; i<size; i++)
        {
            if(!(*(leftAlleles[i]) == *(rightAlleles[i])))
                return false;
        }
        return true;
    }
};


#endif //MUTECT2CPP_MASTER_LOCATIONANDALLELES_H
