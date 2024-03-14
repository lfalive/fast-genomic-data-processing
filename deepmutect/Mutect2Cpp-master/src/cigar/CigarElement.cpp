//
// Created by 梦想家xixi on 2021/11/9.
//

#include "CigarElement.h"
#include "Mutect2Utils.h"

CigarElement::CigarElement(int length, CigarOperator cigarOperator) : length(length), cigarOperator(cigarOperator){
    Mutect2Utils::validateArg(length >= 0, "Cigar element being constructed with negative length");
}

bool CigarElement::operator<(const CigarElement &other) const {
    if(length < other.length)
        return true;
    if(length > other.length)
        return false;
    if(cigarOperator > other.cigarOperator)
        return true;
    if(cigarOperator < other.cigarOperator)
        return false;
    return true;
}

CigarElement::CigarElement() {}

