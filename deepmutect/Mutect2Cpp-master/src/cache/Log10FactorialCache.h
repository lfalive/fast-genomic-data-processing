//
// Created by 梦想家xixi on 2021/11/1.
//

#ifndef MUTECT2CPP_MASTER_LOG10FACTORIALCACHE_H
#define MUTECT2CPP_MASTER_LOG10FACTORIALCACHE_H

#include "IntToDoubleFunctionCache.h"

class Log10FactorialCache : public IntToDoubleFunctionCache{
private:
    static const int CACHE_SIZE = 10000;

protected:
    int maxSize() override {return CACHE_SIZE;}

    double compute(int n) override;
};


#endif //MUTECT2CPP_MASTER_LOG10FACTORIALCACHE_H
