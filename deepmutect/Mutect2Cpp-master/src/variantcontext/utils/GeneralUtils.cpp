//
// Created by 梦想家xixi on 2021/11/27.
//

#include "GeneralUtils.h"
#include <cmath>

double *GeneralUtils::normalizeFromLog10(double *array, int length, bool takeLog10OfOutput, bool keepInLogSpace) {
    double maxValue = arrayMax(array, length);
    if(keepInLogSpace) {
        for(int i = 0; i < length; ++i) {
            array[i] -= maxValue;
        }
        return array;
    } else {
        double * normalized = new double[length];
        for(int i = 0; i < length; ++i) {
            normalized[i] = std::pow(10.0, array[i] - maxValue);
        }
        double sum = 0.0;
        int i;
        for(i = 0; i < length; ++i) {
            sum += normalized[i];
        }
        for(i = 0; i < length; ++i) {
            double x = normalized[i] / sum;
            if (takeLog10OfOutput) {
                x = std::log10(x);
                if (x < -1000000.0 || std::isinf(x)) {
                    x = array[i] - maxValue;
                }
            }

            normalized[i] = x;
        }

        return normalized;
    }
}

double GeneralUtils::arrayMax(double *array, int length) {
    double ret = array[0];
    for(int i = 1; i < length; i++) {
        if(array[i] > ret)
            ret = array[i];
    }
    return ret;
}

double *GeneralUtils::normalizeFromLog10(double *array, int length) {
    return normalizeFromLog10(array, length, false);
}

double *GeneralUtils::normalizeFromLog10(double *array, int length, bool takeLog10OfOutput) {
    return normalizeFromLog10(array, length, takeLog10OfOutput, false);
}
