//
// Created by 梦想家xixi on 2021/10/30.
//

#include "QualityUtils.h"
#include "Mutect2Utils.h"
#include <cmath>
#include <assert.h>

double QualityUtils::errorProbabilityByPhredScore[101]{0};
double QualityUtils::qualToErrorProbCache[255] {0};

uint8_t QualityUtils::errorProbToQual(double errorRate, uint8_t maxQual)
{
    assert(errorRate >= 0.0 && errorRate <= 1.0);
    double d = round(-10.0 * log10(errorRate));
    int qual = std::isinf(d) ? INT32_MAX : (int)d;   // if d is infinity, (int)d will be -2147483647

    return boundQual(qual, maxQual);
}

uint8_t QualityUtils::boundQual(int qual, uint8_t maxQual) {
    return (uint8_t)(std::max(std::min(qual, maxQual & 0xff), 1) & 0xff);
}

uint8_t QualityUtils::errorProbToQual(double errorRate)
{
    return errorProbToQual(errorRate, MAX_SAM_QUAL_SCORE);
}

void QualityUtils::initial() {
    for(int i = 0; i < 255; i++) {
        qualToErrorProbCache[i] = qualToErrorProb((double) i);
    }

    for(int i=0; i<101; i++)
    {
        errorProbabilityByPhredScore[i] = 1.0 / pow(10.0, i/10.0);
    }
}

double QualityUtils::qualToErrorProb(double qual) {
    Mutect2Utils::validateArg(qual >= 0.0, "Qual must be >= 0.0");
    return std::pow(10.0, qual / -10.0);
}

double QualityUtils::qualToErrorProb(uint8_t qual) {
    return qualToErrorProbCache[(int)qual & 0xff];
}

double QualityUtils::qualToErrorProbLog10(double qual)
{
    assert(qual >= 0.0);
    return qual / (-10.0);
}

int QualityUtils::getPhredScoreFromErrorProbability(double probability)
{
    return (int) round(-10.0 * log10(probability));
}