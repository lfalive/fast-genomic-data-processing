//
// Created by 梦想家xixi on 2021/11/12.
//

#include "SWNativeAlignerWrapper.h"
#include "SWNativeResultWrapper.h"
#include "Mutect2Utils.h"
#include "IntelSmithWaterman.h"

SmithWatermanAlignment *SWNativeAlignerWrapper::align(std::shared_ptr<uint8_t[]> reference, int refLength, std::shared_ptr<uint8_t[]> alternate, int altLength, SWParameters *parameters,
                                                      SWOverhangStrategy overhangStrategy) {
    Mutect2Utils::validateArg(parameters, "Null is not allowed there");
    //Mutect2Utils::validateArg(overhangStrategy, "Null is not allowed there");

    int matchIndex = -1;
    if(overhangStrategy == SWOverhangStrategy::SOFTCLIP || overhangStrategy == SWOverhangStrategy::IGNORE) {
        matchIndex = Mutect2Utils::lastIndexOf(reference, refLength, alternate, altLength);
    }

    SmithWatermanAlignment* alignmentResult;
    if(matchIndex != -1) {
        std::vector<CigarElement> lce;
        lce.emplace_back(CigarElement(altLength, M));
        alignmentResult = new SWNativeResultWrapper(new Cigar(lce), matchIndex);
    } else {
        SWNativeAlignerResult* alignment = IntelSmithWaterman::align(reference.get(), refLength, alternate.get(), altLength, parameters,overhangStrategy);
        alignmentResult = new SWNativeResultWrapper(*alignment);
        delete alignment;
    }
    return alignmentResult;
}
