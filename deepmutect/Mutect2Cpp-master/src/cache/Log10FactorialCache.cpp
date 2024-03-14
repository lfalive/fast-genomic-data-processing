//
// Created by 梦想家xixi on 2021/11/1.
//

#include "Log10FactorialCache.h"
#include "DigammaCache.h"
#include <cmath>

double Log10FactorialCache::compute(int n) {
    return DigammaCache::logGamma(n + 1) * 0.4342944819032518;
}
