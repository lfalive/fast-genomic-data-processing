//
// Created by 梦想家xixi on 2021/10/12.
//
#include <iostream>
#include <stdexcept>
#include <cmath>

#include "Mutect2Utils.h"
#include "MathUtils.h"
#include "QualityUtils.h"
#include "NaturalLogUtils.h"

std::string Mutect2Utils::replaceWith(std::string& str1, const std::string& str2, const std::string& str3){
    int pos;
    pos = str1.find(str2);
    while(pos != -1){
        str1.replace(pos,str1.length(),str3);
        pos = str1.find(str2);
    }
    std::cout << str1 << std::endl;
    return str1;
}

bool Mutect2Utils::overlaps(int start, int end, int start2, int end2) {
    return start2 >= start && start2 <= end || end2 >= start && end2 <= end || encloses(start2, end2, start, end);
}

bool Mutect2Utils::encloses(int outerStart, int outerEnd, int innerStart, int innerEnd) {
    return innerStart >= outerStart && innerEnd <= outerEnd;
}

void Mutect2Utils::validateArg(bool condition, const std::string& msg) {
    if(!condition){
        throw std::invalid_argument(msg);
    }
}

bool Mutect2Utils::goodProbability(const double result) {
    return result >= 0.0 && result <= 1.0;
}

double Mutect2Utils::logLikelihoodRatio(const int nRef, const std::vector<uint8_t>& altQuals, const int repeatFactor) {
    int nAlt = repeatFactor * altQuals.size();
    int n = nRef + nAlt;

    double fTildeRatio = std::exp(MathUtils::digamma(nRef + 1) - MathUtils::digamma(nAlt + 1));
    double betaEntropy = MathUtils::log10ToLog(-MathUtils::log10Factorial(n+1) + MathUtils::log10Factorial(nAlt) + MathUtils::log10Factorial(nRef));

    double readSum = 0;
    for(uint8_t qual : altQuals) {
        double epsilon = QualityUtils::qualToErrorProb(qual);
        double zBarAlt = (1 - epsilon) / (1 - epsilon + epsilon * fTildeRatio);
        double logEpsilon = NaturalLogUtils::qualToLogErrorProb(qual);
        double logOneMinusEpsilon = NaturalLogUtils::qualToLogProb(qual);
        readSum += zBarAlt * (logOneMinusEpsilon - logEpsilon) + MathUtils::fastBernoulliEntropy(zBarAlt);
    }

    return betaEntropy + readSum * repeatFactor;

}

double Mutect2Utils::logLikelihoodRatio(int refCount, int altCount, double errorProbability) {
    uint8_t qual = QualityUtils::errorProbToQual(errorProbability);
    std::vector<uint8_t> tmp;
    tmp.emplace_back(qual);
    return logLikelihoodRatio(refCount, tmp, altCount);
}

int Mutect2Utils::lastIndexOf(const std::shared_ptr<uint8_t[]>& reference_, int refLength, const std::shared_ptr<uint8_t[]>& query_, int queryLength) {
    uint8_t * reference = reference_.get();
    uint8_t * query = query_.get();
    for (int r = refLength - queryLength; r >= 0; r--) {
        int q = 0;
        while (q < queryLength && reference[r+q] == query[q]) {
            q++;
        }
        if (q == queryLength) {
            return r;
        }
    }
    return -1;
}

//需要delete
std::shared_ptr<uint8_t[]> Mutect2Utils::copyOfRange(const std::shared_ptr<uint8_t[]>&original, int originalLength, int from, int to, int &length) {
    int newLength = to - from;
    if(newLength < 0)
        throw std::invalid_argument("from > to");
    std::shared_ptr<uint8_t[]> copy(new uint8_t[newLength]);
    length = newLength;
    memcpy(copy.get(), original.get()+from, std::min(originalLength - from, newLength));
    return copy;
}

int Mutect2Utils::Int_compare(int x, int y) {
    return (x < y) ? -1 : ((x == y) ? 0 : 1);
}

uint8_t Mutect2Utils::decodeBase(uint8_t i) {
    switch (i) {
        case 1: {
            return 'A';
        }
        case 2: {
            return 'C';
        }
        case 4: {
            return 'G';
        }
        case 8: {
            return 'T';
        }
        default: {
            return 'N';
        }
    }
}








