//
// Created by 梦想家xixi on 2021/11/13.
//

#include "SWNativeResultWrapper.h"
#include "cigar/TextCigarCodec.h"

SWNativeResultWrapper::SWNativeResultWrapper(const SWNativeAlignerResult& nativeResult) : cigar(TextCigarCodec::decode(nativeResult.cigar)), alignmentOffset(nativeResult.alignment_offset){
}

SWNativeResultWrapper::SWNativeResultWrapper(Cigar *cigar, int alignmentOffset) : cigar(cigar), alignmentOffset(alignmentOffset){}

std::shared_ptr<Cigar> & SWNativeResultWrapper::getCigar() {
    return cigar;
}

int SWNativeResultWrapper::getAlignmentOffset() {
    return alignmentOffset;
}

