//
// Created by lhh on 5/23/22.
//

#include "LocationAndAlleles.h"

LocationAndAlleles::LocationAndAlleles(int loc, std::vector<std::shared_ptr<Allele>> & alleles): loc(loc), alleles(alleles)
{
}

int LocationAndAlleles::getLoc() {
    return loc;
}

std::vector<std::shared_ptr<Allele>> &LocationAndAlleles::getAlleles() {
    return alleles;
}