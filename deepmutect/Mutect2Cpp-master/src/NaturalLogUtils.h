//
// Created by 梦想家xixi on 2021/11/1.
//

#ifndef MUTECT2CPP_MASTER_NATURALLOGUTILS_H
#define MUTECT2CPP_MASTER_NATURALLOGUTILS_H

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

class NaturalLogUtils {
private:

    static constexpr double PHRED_TO_LOG_ERROR_PROB_FACTOR = -0.23025850929940458;
    static constexpr double LOG1MEXP_THRESHOLD = -0.6931471805599453;
    static double qualToLogProbCache[255];

public:
    static double LOG_ONE_HALF;

    // initialization of qualToLogProbCache[], you need to call this method before using any method of this class
    static void initial();
    static double qualToLogErrorProb(double qual);
    static double qualToLogErrorProb(uint8_t qual);
    static double log1mexp(double a);
    static double qualToLogProb(uint8_t qual);

    static std::shared_ptr<std::vector<double>> normalizeLog(std::shared_ptr<std::vector<double>> array, bool takeLogOfOutput, bool inPlace);
    static std::vector<double> normalizeLog(std::vector<double> array, bool takeLogOfOutput, bool inPlace);

    /**
     * Computes $\log(\sum_i e^{a_i})$ trying to avoid underflow issues by using the log-sum-exp trick.
     *
     * <p>
     * This trick consists of shifting all the log values by the maximum so that exponent values are
     * much larger (close to 1) before they are summed. Then the result is shifted back down by
     * the same amount in order to obtain the correct value.
     * </p>
     * @return any double value.
     */
    static double logSumExp(std::shared_ptr<std::vector<double>> logValues);

    static double logSumExp(std::vector<double>& logValues);

    /**
     * normalizes the log-probability array in-place.
     *
     * @param array the array to be normalized
     * @return the normalized-in-place array
     */
    static std::shared_ptr<std::vector<double>> normalizeFromLogToLinearSpace(std::shared_ptr<std::vector<double>> array);

    static std::shared_ptr<std::vector<double>> posteriors( const std::vector<double>& logPriors, const std::vector<double>& logLikelihoods);

    static std::shared_ptr<std::vector<double>> posteriors(std::shared_ptr<std::vector<double>> logPriors, std::shared_ptr<std::vector<double>> logLikelihoods);
};


#endif //MUTECT2CPP_MASTER_NATURALLOGUTILS_H
