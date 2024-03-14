//
// Created by 梦想家xixi on 2021/11/9.
//

#ifndef MUTECT2CPP_MASTER_CIGARELEMENT_H
#define MUTECT2CPP_MASTER_CIGARELEMENT_H


#include "CigarOperator.h"

class CigarElement {
private:
    int length{};
    CigarOperator cigarOperator;

public:
    CigarElement(int length, CigarOperator cigarOperator);
    int getLength() const {return length;}
    CigarOperator getOperator() const {return cigarOperator;}

    bool operator<(const CigarElement & other) const;

    CigarElement();
};


#endif //MUTECT2CPP_MASTER_CIGARELEMENT_H
