//
// Created by 梦想家xixi on 2021/11/13.
//

#ifndef MUTECT2CPP_MASTER_INTELSMITHWATERMAN_H
#define MUTECT2CPP_MASTER_INTELSMITHWATERMAN_H


#include "SWNativeAlignerResult.h"
#include "SWParameters.h"
#include "SmithWatermanAligner.h"

class IntelSmithWaterman {
public:
    static SWNativeAlignerResult* align(uint8_t* refArray, int refLength, uint8_t* altArray, int altLength, SWParameters* parameters, SWOverhangStrategy overhangStrategy);
    static int getStrategy(SWOverhangStrategy strategy);
    static std::string trim(std::string &s);
};


#endif //MUTECT2CPP_MASTER_INTELSMITHWATERMAN_H
