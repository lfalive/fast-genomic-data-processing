//
// Created by 梦想家xixi on 2022/1/4.
//

#ifndef MUTECT2CPP_MASTER_REFERENCECONTEXT_H
#define MUTECT2CPP_MASTER_REFERENCECONTEXT_H


#include "SimpleInterval.h"

class ReferenceContext {
private:
    char refBase;
    std::shared_ptr<SimpleInterval> interval;
    std::shared_ptr<uint8_t[]> refCache;
    int len;

public:
    ReferenceContext(std::shared_ptr<SimpleInterval>  interval, char refBase);
    ReferenceContext(const ReferenceContext& other);
    uint8_t getBase();
    const std::shared_ptr<SimpleInterval> & getInterval();
    void setCache(const std::shared_ptr<uint8_t[]> & refCache, int len);
    std::shared_ptr<uint8_t[]> getCache(int & len);
};


#endif //MUTECT2CPP_MASTER_REFERENCECONTEXT_H
