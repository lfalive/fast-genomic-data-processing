//
// Created by cluster on 22-11-14.
//

#include <climits>
#include <cmath>
#include "BinomialDistribution.h"
#include <stdexcept>
#include "DigammaCache.h"
#include "BetaBinomialDistribution.h"

const std::vector<double > BinomialDistribution::EXACT_STIRLING_ERRORS = {0.0, /* 0.0 */
                                                                          0.1534264097200273452913848, /* 0.5 */
                                                                          0.0810614667953272582196702, /* 1.0 */
                                                                          0.0548141210519176538961390, /* 1.5 */
                                                                          0.0413406959554092940938221, /* 2.0 */
                                                                          0.03316287351993628748511048, /* 2.5 */
                                                                          0.02767792568499833914878929, /* 3.0 */
                                                                          0.02374616365629749597132920, /* 3.5 */
                                                                          0.02079067210376509311152277, /* 4.0 */
                                                                          0.01848845053267318523077934, /* 4.5 */
                                                                          0.01664469118982119216319487, /* 5.0 */
                                                                          0.01513497322191737887351255, /* 5.5 */
                                                                          0.01387612882307074799874573, /* 6.0 */
                                                                          0.01281046524292022692424986, /* 6.5 */
                                                                          0.01189670994589177009505572, /* 7.0 */
                                                                          0.01110455975820691732662991, /* 7.5 */
                                                                          0.010411265261972096497478567, /* 8.0 */
                                                                          0.009799416126158803298389475, /* 8.5 */
                                                                          0.009255462182712732917728637, /* 9.0 */
                                                                          0.008768700134139385462952823, /* 9.5 */
                                                                          0.008330563433362871256469318, /* 10.0 */
                                                                          0.007934114564314020547248100, /* 10.5 */
                                                                          0.007573675487951840794972024, /* 11.0 */
                                                                          0.007244554301320383179543912, /* 11.5 */
                                                                          0.006942840107209529865664152, /* 12.0 */
                                                                          0.006665247032707682442354394, /* 12.5 */
                                                                          0.006408994188004207068439631, /* 13.0 */
                                                                          0.006171712263039457647532867, /* 13.5 */
                                                                          0.005951370112758847735624416, /* 14.0 */
                                                                          0.005746216513010115682023589, /* 14.5 */
                                                                          0.005554733551962801371038690 /* 15.0 */};

BinomialDistribution::BinomialDistribution(int trails, double p) : probabilityOfSuccess(p), numberOfTrials(trails){

}

double BinomialDistribution::cumulativeProbability(int x) {
    double ret;
    if (x < 0) {
        ret = 0.0;
    } else if (x >= numberOfTrials) {
        ret = 1.0;
    } else {
        ret = 1.0 - regularizedBeta(probabilityOfSuccess,
                                         x + 1.0, numberOfTrials - x);
    }
    return ret;
}

double BinomialDistribution::regularizedBeta(double x, double a, double b) {
    return regularizedBeta(x, a, b, DEFAULT_EPSILON, INT_MAX);
}

double BinomialDistribution::regularizedBeta(double x, double a, double b, double epsilon, int maxIterations) {
    double ret;

    if (std::isinf(x) ||
            std::isinf(a) ||
            std::isinf(b) ||
        x < 0 ||
        x > 1 ||
        a <= 0 ||
        b <= 0) {
        ret = nan("");
    } else if (x > (a + 1) / (2 + b + a) &&
               1 - x <= (b + 1) / (2 + b + a)) {
        ret = 1 - regularizedBeta(1 - x, b, a, epsilon, maxIterations);
    } else {
        ret = std::exp(a * std::log(x) + b * std::log1p(-x) - std::log(a) - BetaBinomialDistribution::logBeta(a, b)) * 1.0 / evaluate(x, epsilon, maxIterations, a, b);
    }
    return ret;
}


double BinomialDistribution::evaluate(double x, double epsilon, int maxIterations, double _a, double _b) {
    double small = 1e-50;
    double hPrev = getA(0, x);

    // use the value of small as epsilon criteria for zero checks
    if (std::abs(hPrev) < small) {
        hPrev = small;
    }

    int n = 1;
    double dPrev = 0.0;
    double cPrev = hPrev;
    double hN = hPrev;

    while (n < maxIterations) {
        double a = getA(n, x);
        double b = getB(n, x, _a, _b);

        double dN = a + b * dPrev;
        if (std::abs(dN) < small) {
            dN = small;
        }
        double cN = a + b / cPrev;
        if (std::abs(cN) < small) {
            cN = small;
        }

        dN = 1 / dN;
        double deltaN = cN * dN;
        hN = hPrev * deltaN;

        if (std::isinf(hN)) {
            throw std::invalid_argument("");
        }
        if (std::isinf(hN)) {
            throw std::invalid_argument("");
        }
        if (std::abs(deltaN - 1.0) < epsilon) {
            break;
        }

        dPrev = dN;
        cPrev = cN;
        hPrev = hN;
        n++;
    }

    if (n > maxIterations) {
        throw std::invalid_argument("");
    }

    return hN;
}

double BinomialDistribution::getB(int n, double x, double a, double b) {
    {
        double ret;
        double m;
        if (n % 2 == 0) { // even
            m = n / 2.0;
            ret = (m * (b - m) * x) /
                  ((a + (2 * m) - 1) * (a + (2 * m)));
        } else {
            m = (n - 1.0) / 2.0;
            ret = -((a + m) * (a + b + m) * x) /
                  ((a + (2 * m)) * (a + (2 * m) + 1.0));
        }
        return ret;
    }
}

double BinomialDistribution::logProbability(int x) {
    if (numberOfTrials == 0) {
        return (x == 0) ? 0. : (-1.0 / 0.0);
    }
    double ret;
    if (x < 0 || x > numberOfTrials) {
        ret = (-1.0 / 0.0);
    } else {
        ret = logBinomialProbability(x,numberOfTrials, probabilityOfSuccess,1.0 - probabilityOfSuccess);
    }
    return ret;
}

double BinomialDistribution::logBinomialProbability(int x, int n, double p, double q) {
    double ret;
    if (x == 0) {
        if (p < 0.1) {
            ret = -getDeviancePart(n, n * q) - n * p;
        } else {
            ret = n * std::log(q);
        }
    } else if (x == n) {
        if (q < 0.1) {
            ret = -getDeviancePart(n, n * p) - n * q;
        } else {
            ret = n * std::log(p);
        }
    } else {
        ret = getStirlingError(n) - getStirlingError(x) -
              getStirlingError(n - x) - getDeviancePart(x, n * p) -
              getDeviancePart(n - x, n * q);
        double f = ((2*3.1415926) * x * (n - x)) / n;
        ret = -0.5 * std::log(f) + ret;
    }
    return ret;

}

double BinomialDistribution::getDeviancePart(double x, double mu) {
    double ret;
    if (std::abs(x - mu) < 0.1 * (x + mu)) {
        double d = x - mu;
        double v = d / (x + mu);
        double s1 = v * d;
        double s = (1.0 / 0.0);
        double ej = 2.0 * x * v;
        v *= v;
        int j = 1;
        while (s1 != s) {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1);
            ++j;
        }
        ret = s1;
    } else {
        ret = x * std::log(x / mu) + mu - x;
    }
    return ret;
}

double BinomialDistribution::getStirlingError(double z){
    double ret;
    if (z < 15.0) {
        double z2 = 2.0 * z;
        if (std::floor(z2) == z2) {
            ret = EXACT_STIRLING_ERRORS[(int) z2];
        } else {
            ret = DigammaCache::logGamma(z + 1.0) - (z + 0.5) * std::log(z) +
                  z - 0.9189385332046727;
        }
    } else {
        double z2 = z * z;
        ret = (0.083333333333333333333 -
               (0.00277777777777777777778 -
                (0.00079365079365079365079365 -
                 (0.000595238095238095238095238 -
                  0.0008417508417508417508417508 /
                  z2) / z2) / z2) / z2) / z;
    }
    return ret;
}

double BinomialDistribution::probability(int i) {
    double l = logProbability(i);
    return std::isinf(l) ? 0.0 : std::exp(l);
}

