//
// Created by 梦想家xixi on 2021/11/12.
//

#include "TextCigarCodec.h"

std::string TextCigarCodec::encode(Cigar cigar) {
    if(cigar.isEmpty())
        return "*";
    else {
        std::string ret;
        std::vector<CigarElement> cigarList = cigar.getCigarElements();
        for(CigarElement cigarElement : cigarList) {
            int length = cigarElement.getLength();
            ret += std::to_string(length);
            switch (cigarElement.getOperator()) {
                case M : ret += 'M';
                case I : ret += 'I';
                case D : ret += 'D';
                case N : ret += 'N';
                case S : ret += 'S';
                case H : ret += 'H';
                case P : ret += 'P';
                case X : ret += 'X';
                case EQ : ret += '=';
                default:
                    continue;
            }
        }
        return ret;
    }
}

Cigar *TextCigarCodec::decode(std::string textCigar) {
    if("*" == textCigar)
        return new Cigar();
    else {
        Cigar* ret = new Cigar();
        for(int i = 0; i < textCigar.size(); i++) {
            if(!isDigit(textCigar[i])) {
                throw std::invalid_argument("Malformed CIGAR string");
            }
            int length = textCigar[i] - ZERO_BYTE;
            ++i;

            while(isDigit(textCigar[i])) {
                length = length * 10 + textCigar[i] - ZERO_BYTE;
                ++i;
            }
            CigarOperator cigarOperator = CigarOperatorUtils::characterToEnum(textCigar[i]);
            ret->add(CigarElement(length, cigarOperator));
        }
        return ret;
    }
}

bool TextCigarCodec::isDigit(uint8_t c) {
    return c >= ZERO_BYTE && c <= NINE_BYTE;
}
