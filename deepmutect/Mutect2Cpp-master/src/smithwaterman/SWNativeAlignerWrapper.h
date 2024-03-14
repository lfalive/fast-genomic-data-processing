//
// Created by 梦想家xixi on 2021/11/12.
//

#ifndef MUTECT2CPP_MASTER_SWNATIVEALIGNERWRAPPER_H
#define MUTECT2CPP_MASTER_SWNATIVEALIGNERWRAPPER_H

#include "SmithWatermanAligner.h"
#include "SmithWatermanAlignment.h"
#include "SWNativeAlignerResult.h"
#include "cigar/TextCigarCodec.h"

class SWNativeAlignerWrapper : public SmithWatermanAligner{
public:
    SmithWatermanAlignment* align(std::shared_ptr<uint8_t[]> ref, int refLength, std::shared_ptr<uint8_t[]> alt, int altLength, SWParameters* parameters, SWOverhangStrategy overhangStrategy);
};


#endif //MUTECT2CPP_MASTER_SWNATIVEALIGNERWRAPPER_H
