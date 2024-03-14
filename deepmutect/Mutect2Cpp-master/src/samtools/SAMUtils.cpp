//
// Created by 梦想家xixi on 2021/12/20.
//

#include <iostream>
#include "SAMUtils.h"

int SAMUtils::getUnclippedStart(int alignmentStart, std::shared_ptr<Cigar> cigar) {
    int unClippedStart = alignmentStart;
    for(CigarElement cig : cigar->getCigarElements()) {
        CigarOperator op = cig.getOperator();
        if (op != S && op != H) {
            break;
        }
        unClippedStart -= cig.getLength();
    }
    return unClippedStart;
}

int SAMUtils::getUnclippedEnd(int alignmentEnd, std::shared_ptr<Cigar> cigar) {
    int unClippedEnd = alignmentEnd;
    std::vector<CigarElement> cigs = cigar->getCigarElements();
    for(int i = cigs.size() - 1; i >= 0; --i) {
        CigarElement cig = cigs[i];
        CigarOperator op = cig.getOperator();
        if (op != S && op != H) {
            break;
        }

        unClippedEnd += cig.getLength();
    }

    return unClippedEnd;

}

bool SAMUtils::isValidUnsignedIntegerAttribute(long value) {
    return value >= 0L && value <= 4294967295L;
}

short SAMUtils::makeBinaryTag(std::string &tag) {
    if(tag.length() != 2) {
        throw std::invalid_argument("String tag does not have length() == 2: ");
    }
    char a = tag[1];
    char b = tag[0];
    return a << 8 | b;
}

std::shared_ptr<uint8_t[]> SAMUtils::fastqToPhred(std::string &fastq) {
    if(fastq.empty()) {
        return nullptr;
    } else {
        int length = fastq.size();
        std::shared_ptr<uint8_t[]> scores(new uint8_t[length]);
        uint8_t * scores_ = scores.get();
        for(int i = 0; i < length; ++i) {
            scores_[i] = (uint8_t)fastqToPhred(fastq[i]);
        }
        return scores;
    }
}

int SAMUtils::fastqToPhred(char ch) {
    if(ch >= '!' && ch <= '~') {
        return ch - 33;
    } else {
        throw std::invalid_argument("Invalid fastq character");
    }
}

char SAMUtils::phredToFastq(int phredScore) {
    if(phredScore >= 0 && phredScore <= 93) {
        return (char)(33 + phredScore);
    } else {
        throw std::invalid_argument("Cannot encode phred score");
    }
}

std::string SAMUtils::phredToFastq(std::shared_ptr<uint8_t[]>buffer, int offset, int length) {
    char* chars = new char[length];
    uint8_t * buffer_ = buffer.get();
    for(int i = 0; i < length; ++i) {
        chars[i] = phredToFastq(buffer_[offset + i] & 255);
    }

    std::string ret(chars, length);
    delete[] chars;
    return ret;
}

std::string SAMUtils::phredToFastq(std::shared_ptr<uint8_t[]>buffer, int length) {
    if(buffer == nullptr)
        return "";
    else
        return phredToFastq(buffer, 0, length);
}
