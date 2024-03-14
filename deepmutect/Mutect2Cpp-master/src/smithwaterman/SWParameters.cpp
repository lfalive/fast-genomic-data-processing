//
// Created by 梦想家xixi on 2021/11/12.
//

#include "SWParameters.h"
#include "Mutect2Utils.h"

SWParameters::SWParameters(int matchValue, int mismatchPenalty, int gapOpenPenalty, int gapExtendPenalty) : matchValue(matchValue), mismatchPenalty(mismatchPenalty), gapOpenPenalty(gapOpenPenalty), gapExtendPenalty(gapExtendPenalty){
    Mutect2Utils::validateArg(matchValue >= 0, "matchValue must be >= 0");
    Mutect2Utils::validateArg(mismatchPenalty <= 0, "mismatchPenalty must be <= 0");
    Mutect2Utils::validateArg(gapOpenPenalty <= 0, "gapOpenPenalty must be <= 0");
    Mutect2Utils::validateArg(gapExtendPenalty <= 0, "gapExtendPenalty must be <= 0");
}
