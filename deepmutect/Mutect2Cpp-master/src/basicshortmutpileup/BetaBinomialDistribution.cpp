//
// Created by cluster on 22-11-8.
//

#include "BetaBinomialDistribution.h"
#include "CombinatoricsUtils.h"
#include "boost/math/distributions.hpp"
#include <stdexcept>
#include <cfloat>
#include "cache/DigammaCache.h"

std::vector<double> BetaBinomialDistribution::DELTA = {
        .833333333333333333333333333333E-01,
        -.277777777777777777777777752282E-04,
        .793650793650793650791732130419E-07,
        -.595238095238095232389839236182E-09,
        .841750841750832853294451671990E-11,
        -.191752691751854612334149171243E-12,
        .641025640510325475730918472625E-14,
        -.295506514125338232839867823991E-15,
        .179643716359402238723287696452E-16,
        -.139228964661627791231203060395E-17,
        .133802855014020915603275339093E-18,
        -.154246009867966094273710216533E-19,
        .197701992980957427278370133333E-20,
        -.234065664793997056856992426667E-21,
        .171348014966398575409015466667E-22
};

BetaBinomialDistribution::BetaBinomialDistribution(double alpha, double beta, int n, int rng) : alpha(alpha), beta(beta), n(n), rng(rng){
    if(alpha < 0 || beta < 0 || n < 0) {
        std::string param = ("alpha, beta and n must be greater than zero.");
        throw std::invalid_argument(param);
    }
}

double BetaBinomialDistribution::logProbability(int k) {
    if(k < 0) {
        std::string param = ("alpha, beta and n must be greater than zero.");
        throw std::invalid_argument(param);
    }
    double a = logBeta(k+alpha, n - k + beta);
    double b = logBeta(alpha, beta);
    double c = logBeta(1, 2);
    return k > n ? -DBL_MAX : CombinatoricsUtils::binomialCoefficientLog(n, k) + logBeta(k+alpha, n - k + beta) - logBeta(alpha, beta);
}

double BetaBinomialDistribution::logBeta(double p, double q) {
    double a = std::min(p, q);
    double b = std::max(p, q);
    if (a >= 10.0) {
        double w = sumDeltaMinusDeltaSum(a, b);
        double h = a / b;
        double c = h / (1.0 + h);
        double u = -(a - 0.5) * std::log(c);
        double v = b * std::log1p(h);
        if (u <= v) {
            return (((-0.5 * std::log(b) + HALF_LOG_TWO_PI) + w) - u) - v;
        } else {
            return (((-0.5 * std::log(b) + HALF_LOG_TWO_PI) + w) - v) - u;
        }
    } else if (a > 2.0) {
        if (b > 1000.0) {
            int n = (int) std::floor(a - 1.0);
            double prod = 1.0;
            double ared = a;
            for (int i = 0; i < n; i++) {
                ared -= 1.0;
                prod *= ared / (1.0 + ared / b);
            }
            return (std::log(prod) - n * std::log(b)) +
                   (DigammaCache::logGamma(ared) +
                    logGammaMinusLogGammaSum(ared, b));
        } else {
            double prod1 = 1.0;
            double ared = a;
            while (ared > 2.0) {
                ared -= 1.0;
                double h = ared / b;
                prod1 *= h / (1.0 + h);
            }
            if (b < 10.0) {
                double prod2 = 1.0;
                double bred = b;
                while (bred > 2.0) {
                    bred -= 1.0;
                    prod2 *= bred / (ared + bred);
                }
                return std::log(prod1) +
                        std::log(prod2) +
                       (DigammaCache::logGamma(ared) +
                        (DigammaCache::logGamma(bred) -
                         logGammaSum(ared, bred)));
            } else {
                return std::log(prod1) +
                        DigammaCache::logGamma(ared) +
                       logGammaMinusLogGammaSum(ared, b);
            }
        }
    } else if (a >= 1.0) {
        if (b > 2.0) {
            if (b < 10.0) {
                double prod = 1.0;
                double bred = b;
                while (bred > 2.0) {
                    bred -= 1.0;
                    prod *= bred / (a + bred);
                }
                return std::log(prod) +
                       (DigammaCache::logGamma(a) +
                        (DigammaCache::logGamma(bred) -
                         logGammaSum(a, bred)));
            } else {
                return DigammaCache::logGamma(a) +
                       logGammaMinusLogGammaSum(a, b);
            }
        } else {
            return DigammaCache::logGamma(a) +
                    DigammaCache::logGamma(b) -
                   logGammaSum(a, b);
        }
    } else {
        if (b >= 10.0) {
            return DigammaCache::logGamma(a) +
                   logGammaMinusLogGammaSum(a, b);
        } else {
            // The following command is the original NSWC implementation.
            // return Gamma.logGamma(a) +
            // (Gamma.logGamma(b) - Gamma.logGamma(a + b));
            // The following command turns out to be more accurate.
            return std::log(DigammaCache::gamma(a) * DigammaCache::gamma(b) /
                                    DigammaCache::gamma(a + b));
        }
    }
}

double BetaBinomialDistribution::sumDeltaMinusDeltaSum(double p, double q) {
    double a = std::min(p, q);
    double b = std::max(p, q);
    double sqrtT = 10.0 / a;
    double t = sqrtT * sqrtT;
    double z = DELTA[DELTA.size() - 1];
    for (int i = DELTA.size() - 2; i >= 0; i--) {
        z = t * z + DELTA[i];
    }
    return z / a + deltaMinusDeltaSum(a, b);
}

double BetaBinomialDistribution::deltaMinusDeltaSum(double a, double b) {
    double h = a / b;
    double p = h / (1.0 + h);
    double q = 1.0 / (1.0 + h);
    double q2 = q * q;
    /*
     * s[i] = 1 + q + ... - q**(2 * i)
     */
    std::vector<double> s(DELTA.size(), 0);
    s[0] = 1.0;
    for (int i = 1; i < s.size(); i++) {
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }
    /*
     * w = Delta(b) - Delta(a + b)
     */
    double sqrtT = 10.0 / b;
    double t = sqrtT * sqrtT;
    double w = DELTA[DELTA.size() - 1] * s[s.size() - 1];
    for (int i = DELTA.size() - 2; i >= 0; i--) {
        w = t * w + DELTA[i] * s[i];
    }
    return w * p / b;
}

double BetaBinomialDistribution::logGammaMinusLogGammaSum(double a, double b) {
    double d;
    double w;
    if (a <= b) {
        d = b + (a - 0.5);
        w = deltaMinusDeltaSum(a, b);
    } else {
        d = a + (b - 0.5);
        w = deltaMinusDeltaSum(b, a);
    }

    double u = d * std::log1p(a / b);
    double v = a * (std::log(b) - 1.0);

    return u <= v ? (w - u) - v : (w - v) - u;
}

double BetaBinomialDistribution::logGammaSum(double  a, double b) {
    double x = (a - 1.0) + (b - 1.0);
    if (x <= 0.5) {
        return DigammaCache::logGamma1p(1.0 + x);
    } else if (x <= 1.5) {
        return DigammaCache::logGamma1p(x) + std::log1p(x);
    } else {
        return DigammaCache::logGamma1p(x - 1.0) + std::log(x * (1.0 + x));
    }
}

