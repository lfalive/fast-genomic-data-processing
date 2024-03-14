//
// Created by 梦想家xixi on 2021/11/12.
//

#ifndef MUTECT2CPP_MASTER_TEXTCIGARCODEC_H
#define MUTECT2CPP_MASTER_TEXTCIGARCODEC_H


#include <cstdint>
#include <string>
#include "cigar/Cigar.h"

class TextCigarCodec {
private:
    static const uint8_t ZERO_BYTE = '0';
    static const uint8_t NINE_BYTE = '9';
    static bool isDigit(uint8_t);

public:
    TextCigarCodec() = default;

    static std::string encode(Cigar cigar);
    static Cigar* decode(std::string textCigar);
};


#endif //MUTECT2CPP_MASTER_TEXTCIGARCODEC_H
