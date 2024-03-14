//
// Created by 梦想家xixi on 2022/1/4.
//

#include "ReferenceContext.h"
#include <utility>

ReferenceContext::ReferenceContext(std::shared_ptr<SimpleInterval>  interval, char refBase) : interval(std::move(interval)), refBase(refBase), refCache(nullptr), len(0){

}

ReferenceContext::ReferenceContext(const ReferenceContext &other) : interval(other.interval), refBase(other.refBase), refCache(other.refCache), len(other.len) {
}

uint8_t ReferenceContext::getBase() {
    return refBase;
}

const std::shared_ptr<SimpleInterval> & ReferenceContext::getInterval() {
    return interval;
}

void ReferenceContext::setCache(const std::shared_ptr<uint8_t[]> & refCache, int len) {
    this->refCache = refCache;
    this->len = len;
}

std::shared_ptr<uint8_t[]> ReferenceContext::getCache(int &len) {
    len = this->len;
    return refCache;
}
