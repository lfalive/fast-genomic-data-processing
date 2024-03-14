//
// Created by 梦想家xixi on 2021/10/30.
//

#ifndef MUTECT2CPP_MASTER_QUALITYUTILS_H
#define MUTECT2CPP_MASTER_QUALITYUTILS_H


#include <cstdint>

class QualityUtils {
private:
    static double errorProbabilityByPhredScore[101];
    static double qualToErrorProbCache[255];

public:
    const static char MIN_USABLE_Q_SCORE = 6;
    const static int MAPPING_QUALITY_UNAVALIABLE = 255;

    /**
     * bams containing quals above this value are extremely suspicious and we should warn the user
     */
    const static char MAX_REASONABLE_Q_SCORE = 60;

    const static char MAX_SAM_QUAL_SCORE = 93;

    /**
     * Convert a probability of being wrong to a phred-scaled quality score (0.01 => 20).
     *
     * Note, this function caps the resulting quality score by the public static value MIN_REASONABLE_ERROR
     * and by 1 at the low-end.
     *
     * WARNING -- because this function takes a byte for maxQual, you must be careful in converting
     * integers to byte.  The appropriate way to do this is ((byte)(myInt & 0xFF))
     *
     * @param errorRate a probability (0.0-1.0) of being wrong (i.e., 0.01 is 1% change of being wrong)
     * @return a quality score (0-maxQual)
     */
    static uint8_t errorProbToQual(double prob, uint8_t maxQual);

    /**
      * Convert a probability of being wrong to a phred-scaled quality score (0.01 => 20).
      *
      * Note, this function caps the resulting quality score by the public static value MAX_SAM_QUAL_SCORE
      * and by 1 at the low-end.
      *
      * @param errorRate a probability (0.0-1.0) of being wrong (i.e., 0.01 is 1% change of being wrong)
      * @return a quality score (0-MAX_SAM_QUAL_SCORE)
      */
    static uint8_t errorProbToQual(double errorRate);
    static uint8_t boundQual(int qual, uint8_t maxQual);
    static double qualToErrorProb(double qual);
    static double qualToErrorProb(uint8_t qual);

    /**
     * Convert a phred-scaled quality score to its log10 probability of being wrong (Q30 => log10(0.001))
     *
     * This is the Phred-style conversion, *not* the Illumina-style conversion.
     *
     * The calculation is extremely efficient
     *
     * @param qual a phred-scaled quality score encoded as a double
     * @return log of probability (0.0-1.0)
     */
    static double qualToErrorProbLog10(double qual);

    static void initial();


    /** Gets the phred score for any given probability of error. */
    static int getPhredScoreFromErrorProbability(double probability);


};


#endif //MUTECT2CPP_MASTER_QUALITYUTILS_H
