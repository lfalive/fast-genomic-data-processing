//
// Created by 梦想家xixi on 2021/11/27.
//

#ifndef MUTECT2CPP_MASTER_GENERALUTILS_H
#define MUTECT2CPP_MASTER_GENERALUTILS_H


class GeneralUtils {
public:
    static double* normalizeFromLog10(double * array, int length, bool takeLog10OfOutput, bool keepInLogSpace);
    static double* normalizeFromLog10(double * array, int length);
    static double* normalizeFromLog10(double * array, int length, bool takeLog10OfOutput);
    static double arrayMax(double * array, int length);
};


#endif //MUTECT2CPP_MASTER_GENERALUTILS_H
