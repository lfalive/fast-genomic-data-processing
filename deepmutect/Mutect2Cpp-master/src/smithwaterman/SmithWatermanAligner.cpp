//
// Created by 梦想家xixi on 2021/11/12.
//

#include "SmithWatermanAligner.h"
#include "SWNativeAlignerWrapper.h"

SWParameters SmithWatermanAligner::ORIGINAL_DEFAULT = SWParameters(3, -1, -4, -3);
SWParameters SmithWatermanAligner::STANDARD_NGS = SWParameters(25, -50, -110, -6);

SmithWatermanAligner * SmithWatermanAligner::getAligner(Implementation type)
{
    switch (type) {
        case FASTEST_AVAILABLE:
            return new SWNativeAlignerWrapper;
        default:
            return nullptr;
    }
}