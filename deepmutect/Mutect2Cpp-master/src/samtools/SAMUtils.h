//
// Created by 梦想家xixi on 2021/12/20.
//

#ifndef MUTECT2CPP_MASTER_SAMUTILS_H
#define MUTECT2CPP_MASTER_SAMUTILS_H


#include "Cigar.h"
#include <string>
#include <memory>

class SAMUtils {
public:
    static int getUnclippedStart(int alignmentStart, std::shared_ptr<Cigar> cigar);
    static int getUnclippedEnd(int alignmentEnd, std::shared_ptr<Cigar> cigar);
    static bool isValidUnsignedIntegerAttribute(long value);
    static short makeBinaryTag(std::string& tag);
    static std::shared_ptr<uint8_t[]> fastqToPhred(std::string& fastq);
    static int fastqToPhred(char ch);
    static char phredToFastq(int phredScore);
    static std::string phredToFastq(std::shared_ptr<uint8_t[]> buffer, int offset, int length);
    static std::string phredToFastq(std::shared_ptr<uint8_t[]> buffer, int length);
};


#endif //MUTECT2CPP_MASTER_SAMUTILS_H
