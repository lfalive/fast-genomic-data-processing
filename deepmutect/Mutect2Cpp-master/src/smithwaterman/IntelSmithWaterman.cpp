//
// Created by 梦想家xixi on 2021/11/13.
//

#include "IntelSmithWaterman.h"
#include "intel/smithwaterman/IntelSmithWaterman.h"

SWNativeAlignerResult*
IntelSmithWaterman::align(uint8_t *refArray, int refLength, uint8_t *altArray, int altLength, SWParameters* parameters,
                          SWOverhangStrategy overhangStrategy) {
    int intStrategy = getStrategy(overhangStrategy);
    int cigarLength = 2*std::max(refLength, altLength);
    auto * cigar = new uint8_t[cigarLength];
	std::fill_n(cigar, cigarLength, 0);
    int offset = SmithWaterman_align(refArray, refLength, altArray, altLength, cigar, cigarLength, parameters->getMatchValue(), parameters->getMismatchPenalty(),parameters->getGapOpenPenalty(), parameters->getGapExtendPenalty(), intStrategy);
    std::string ret = std::string((char*)cigar);
    delete[] cigar;
    return new SWNativeAlignerResult(trim(ret), offset);
}

int IntelSmithWaterman::getStrategy(SWOverhangStrategy strategy) {
    int intStrategy = 0;

    switch(strategy){
        case SOFTCLIP: intStrategy = 9;
            break;
        case INDEL: intStrategy = 10;
            break;
        case LEADING_INDEL: intStrategy = 11;
            break;
        case IGNORE: intStrategy = 12;
            break;
    }

    return intStrategy;
}

std::string IntelSmithWaterman::trim(std::string &s) {
    int len = s.length();
    int st = 0;
    while ((st < len) && (s[st] <= ' ')) {
        st++;
    }
    while ((st < len) && (s[len - 1] <= ' ')) {
        len--;
    }
    return ((st > 0) || (len < s.length())) ? s.substr(st, len) : s;
}
