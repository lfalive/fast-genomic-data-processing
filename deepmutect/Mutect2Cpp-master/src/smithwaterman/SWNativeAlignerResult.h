//
// Created by 梦想家xixi on 2021/11/12.
//

#ifndef MUTECT2CPP_MASTER_SWNATIVEALIGNERRESULT_H
#define MUTECT2CPP_MASTER_SWNATIVEALIGNERRESULT_H

#include <string>
#include <utility>

class SWNativeAlignerResult {
public:
    std::string cigar;
    int alignment_offset;

    SWNativeAlignerResult(std::string cigar, int alignment_offset) : cigar(std::move(cigar)), alignment_offset(alignment_offset){}
};


#endif //MUTECT2CPP_MASTER_SWNATIVEALIGNERRESULT_H
