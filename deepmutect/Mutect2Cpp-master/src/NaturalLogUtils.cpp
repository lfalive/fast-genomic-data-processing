//
// Created by 梦想家xixi on 2021/11/1.
//

#include <cassert>
#include "NaturalLogUtils.h"
#include "Mutect2Utils.h"
#include "MathUtils.h"

double NaturalLogUtils::LOG_ONE_HALF = log(0.5);
double NaturalLogUtils::qualToLogProbCache[255] {0};

double NaturalLogUtils::qualToLogErrorProb(double qual) {
    Mutect2Utils::validateArg(qual >= 0.0, "qual must be >= 0.0 ");
    return qual * PHRED_TO_LOG_ERROR_PROB_FACTOR;
}

double NaturalLogUtils::qualToLogErrorProb(uint8_t qual) {
    return qualToLogErrorProb((double)(qual & 0xff));
}

double NaturalLogUtils::log1mexp(double a) {
    if (a > 0) return std::nan("");
    if (a == 0) return -std::numeric_limits<double>::infinity();
    return (a < LOG1MEXP_THRESHOLD) ? std::log1p(-std::exp(a)) : std::log(-std::expm1(a));
}

void NaturalLogUtils::initial() {
    for(int i = 0; i < 255; i++) {
        qualToLogProbCache[i] = log1mexp(qualToLogErrorProb(double(i)));
    }
}

double NaturalLogUtils::qualToLogProb(uint8_t qual) {
    return qualToLogProbCache[(int) qual & 0xff];
}

std::shared_ptr<std::vector<double>>
NaturalLogUtils::normalizeLog(std::shared_ptr<std::vector<double>> array, bool takeLogOfOutput, bool inPlace) {
    double logSum = logSumExp(array);
    std::shared_ptr<std::vector<double>> result = inPlace ? array : make_shared<std::vector<double>>(*array);
    for(int m=0; m<array->size(); m++)
    {
        result->operator[](m) = array->operator[](m) - logSum;
    }
    return takeLogOfOutput ? result : MathUtils::applyToArrayInPlace(result, [](double x){return exp(x);});
}

double NaturalLogUtils::logSumExp(std::shared_ptr<std::vector<double>> logValues) {
    assert(logValues != nullptr);
    return logSumExp(*logValues);
}

double NaturalLogUtils::logSumExp(std::vector<double>& logValues) {
    int maxElementIndex = MathUtils::maxElementIndex(logValues);
    double maxValue = logValues[maxElementIndex];
    if(maxValue == -std::numeric_limits<double>::infinity())
        return maxValue;
    double sum = 1.0;
    for(int i=0; i<logValues.size(); i++)
    {
        double curVal = logValues[i];
        if(i == maxElementIndex || isinf(curVal))
            continue;
        else {
            sum += exp(curVal - maxValue);
        }
    }
    if(isnan(sum) || isinf(sum))
        throw std::invalid_argument("log10 p: Values must be non-infinite and non-NAN");

    return maxValue + (sum != 1.0 ? log(sum) : 0.0);
}

std::shared_ptr<std::vector<double>>
NaturalLogUtils::normalizeFromLogToLinearSpace(std::shared_ptr<std::vector<double>> array) {
    return normalizeLog(array, false, true);
}

std::shared_ptr<std::vector<double>> NaturalLogUtils::posteriors(std::shared_ptr<std::vector<double>> logPriors,
                                                                 std::shared_ptr<std::vector<double>> logLikelihoods) {
    return posteriors(*logPriors, *logLikelihoods);
}

std::shared_ptr<std::vector<double>> NaturalLogUtils::posteriors(const std::vector<double>& logPriors,
                                                                 const std::vector<double>& logLikelihoods) {
    int size = logLikelihoods.size();
    auto addedArray = make_shared<std::vector<double>>(size, 0.0);
    for(int i=0; i<size; i++)
    {
        addedArray->operator[](i) = logPriors[i] + logLikelihoods[i];
    }
    return normalizeFromLogToLinearSpace(addedArray);
}

std::vector<double> NaturalLogUtils::normalizeLog(std::vector<double> array, bool takeLogOfOutput, bool inPlace) {
    double logSum = logSumExp(array);
    std::vector<double> result = inPlace ? array :std::vector<double>(array.size());
    for(int m=0; m<array.size(); m++)
    {
        result[m] = array[m] - logSum;
    }
    if(!takeLogOfOutput) {
       MathUtils::applyToArrayInPlace(result, [](double x){return exp(x);});
    }
    return result;
}


