//
// Created by 梦想家xixi on 2021/11/13.
//

#ifndef MUTECT2CPP_MASTER_SWNATIVERESULTWRAPPER_H
#define MUTECT2CPP_MASTER_SWNATIVERESULTWRAPPER_H

#include "cigar/Cigar.h"
#include "SWNativeAlignerResult.h"
#include "SmithWatermanAlignment.h"

class SWNativeResultWrapper : public SmithWatermanAlignment{
private:
    std::shared_ptr<Cigar> cigar;
    int alignmentOffset;

public:
    SWNativeResultWrapper(const SWNativeAlignerResult& nativeResult);
    SWNativeResultWrapper(Cigar* cigar, int alignmentOffset);
    std::shared_ptr<Cigar> & getCigar();
    int getAlignmentOffset();
};


#endif //MUTECT2CPP_MASTER_SWNATIVERESULTWRAPPER_H
