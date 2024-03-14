//
// Created by 梦想家xixi on 2021/11/19.
//

#include "SplitCommonSuffices.h"
#include "CommonSuffixSplitter.h"

bool SplitCommonSuffices::tryToTransform(std::shared_ptr<SeqVertex> bottom) {
    Mutect2Utils::validateArg(bottom.get(), "null is not allowed there");
    if(alreadySplit.find(bottom) != alreadySplit.end()) {
        return false;
    } else {
        alreadySplit.insert(bottom);
        return CommonSuffixSplitter::split(getGraph(), bottom);
    }
}
