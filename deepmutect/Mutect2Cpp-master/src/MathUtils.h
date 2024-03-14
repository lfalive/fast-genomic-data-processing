/**
 * A helper class to math calculation
 */

#ifndef MUTECT2CPP_MASTER_MATHUTILS_H
#define MUTECT2CPP_MASTER_MATHUTILS_H

#include <vector>
#include <functional>
#include "cache/DigammaCache.h"
#include "cache/Log10FactorialCache.h"
#include "AssemblyRegion.h"
#include "ReadCache.h"


class MathUtils {
private:
    static DigammaCache & DIGAMMA_CACHE();
    static Log10FactorialCache & LOG_10_FACTORIAL_CACHE();
public:
    static double digamma(int n);

    static double log10ToLog(double log10);

    static double log10Factorial(int n);

    static double fastBernoulliEntropy(double p);

    /**
     * Calculate f(x) = Normal(x | mu = mean, sigma = sd)
     * @param mean the desired mean of the Normal distribution
     * @param sd the desired standard deviation of the Normal distribution
     * @param x the value to evaluate
     * @return a well-formed double
     */
    static double normalDistribution(double mean, double sd, double x);

    static std::vector<double> * normalizeSumToZero(std::vector<double> * array);

    // A fast implementation of the Math.round() method.  This method does not perform
    // under/overflow checking, so this shouldn't be used in the general case (but is fine
    // if one is already make those checks before calling in to the rounding).
    static int fastRound(double d);

    /*
     * bionomial Probability(int, int, double ) with log applied to result
     */
    static double log10BinomialProbability(int n, int k, double log10p);

    /*
     *  Calculates the log10 of the gamma function for x
     * @param x
     * @return
     */
    static double log10Gamma(double x);

    static double log10Fractorial(int n);

    static double log10BinomialCoefficient(int n, int k);

    static double sum(vector<double>& values);

    static double sum(shared_ptr<vector<double>> values);

    /**
    * Apply a method for an array of elements
    * Returns a new array -- the original array in not modified.
    *
    */
    static shared_ptr<vector<double>> applyToArray(shared_ptr<vector<double>> array, double (*func)(double));

    static shared_ptr<vector<double>> applyToArray(vector<double>& array, double (*func)(double));

    /**
    * Apply a method for an array of elements
    * the original array in modified in place.
    *
    */
    static shared_ptr<vector<double>> applyToArrayInPlace(shared_ptr<vector<double>> array, function<double(double)> func);

    static void applyToArrayInPlace(vector<double>& array, function<double(double)> func);

    // sum of int -> double[] function mapped to an index range
    static shared_ptr<vector<double>> sumArrayFunction(int min, int max, function<shared_ptr<vector<double>>(int)> func);

    static int maxElementIndex(shared_ptr<vector<double>> array);

    static int maxElementIndex(vector<double>& array);

    static double distance1(vector<double>& p1, vector<double>& p2);

    /**
    * normalizes the real-space probability array.
    *
    * Does not assume anything about the values in the array, beyond that no elements are below 0.  It's ok
    * to have values in the array of > 1, or have the sum go above 0.
    *
    * @param array the array to be normalized
    * @return a newly allocated array corresponding the normalized values in array
    */
    static shared_ptr<vector<double>> normalizeSumToOne(shared_ptr<vector<double>> array);

    static int median(std::vector<int> &values);

    constexpr static double LOG10_OF_E = 0.4342944819032518;

    static std::vector<double> normalizeLog10(std::vector<double>& array, bool normalizeLog10, bool inPlace);

    static double log10SumLog10(const std::vector<double> &log10Values, int start, int finish);

    static double log10SumLog10(const std::vector<double> &log10Values);

    static int maxElementIndex(const std::vector<double> & array, int start, int endIndex);

    static double logToLog10(const double ln);
};


#endif //MUTECT2CPP_MASTER_MATHUTILS_H
