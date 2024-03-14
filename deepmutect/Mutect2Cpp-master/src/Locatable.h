//
// Created by 梦想家xixi on 2021/10/12.
//

#ifndef MUTECT2CPP_MASTER_LOCATABLE_H
#define MUTECT2CPP_MASTER_LOCATABLE_H

#include "Mutect2Utils.h"
#include <string>

class Locatable {
public:
    virtual std::string getContig() const = 0;
    virtual int getStart() const = 0;
    virtual int getEnd() const = 0;

    virtual int getLengthOnReference();
    virtual bool contains(const std::shared_ptr<Locatable> & other);
    virtual bool contigsMatch(const std::shared_ptr<Locatable> & other);
    virtual bool withinDistanceOf(const std::shared_ptr<Locatable> & other, int distance);
    virtual bool overlaps(const std::shared_ptr<Locatable> & other);
};


#endif //MUTECT2CPP_MASTER_LOCATABLE_H
