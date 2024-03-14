//
// Created by 梦想家xixi on 2021/12/6.
//

#ifndef MUTECT2CPP_MASTER_BASEUTILS_H
#define MUTECT2CPP_MASTER_BASEUTILS_H


#include <cstdint>
#include <memory>

class BaseUtils {
private:
    static int baseIndexMap[256];

public:
    static void initial();
    static bool isRegularBase(uint8_t base);
    static int simpleBaseToBaseIndex(uint8_t base);
    static bool isAllRegularBases(const std::shared_ptr<uint8_t[]>& bases, int length);
    static uint8_t getComplement(uint8_t base);
    static uint8_t getUpper(uint8_t base);
};


#endif //MUTECT2CPP_MASTER_BASEUTILS_H
