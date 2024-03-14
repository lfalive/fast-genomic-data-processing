//
// Created by 梦想家xixi on 2021/11/1.
//

#include "DigammaCache.h"
#include <cmath>
#include <stdexcept>

double DigammaCache::LANCZOS[15] = {
0.99999999999999709182,
57.156235665862923517,
-59.597960355475491248,
14.136097974741747174,
-0.49191381609762019978,
.33994649984811888699e-4,
.46523628927048575665e-4,
-.98374475304879564677e-4,
.15808870322491248884e-3,
-.21026444172410488319e-3,
.21743961811521264320e-3,
-.16431810653676389022e-3,
.84418223983852743293e-4,
-.26190838401581408670e-4,
.36899182659531622704e-5,
};

double DigammaCache::compute(int n) {
    return digamma(n);
}

double DigammaCache::digamma(double n) {
    if(n > 0 && n <= 1e-5)
        return -0.577215664901532860606512090082 - 1 / n;

    if(n >= 49.0) {
        double inv = 1.0 / (n * n);
        return std::log(n) - 0.5 / n - inv * ((1.0 / 12) + inv * (1.0 / 120 - inv / 252));
    }

    return digamma(n + 1) - 1 / n;
}

double DigammaCache::logGamma(double x) {
    double ret;
    double HALF_LOG_2_PI = 0.5 * std::log(2.0 * (105414357.0 / 33554432.0 + 1.984187159361080883e-9));
    if(x != x || x <= 0) {
        ret = std::nan("");
    } else if (x < 0.5) {
        return logGamma1p(x) - std::log(x);
    } else if (x <= 2.5) {
        return logGamma1p((x -0.5) - 0.5);
    } else if (x <= 8.0) {
        int n = (int) std::floor(x - 1.5);
        double prod = 1.0;
        for (int i = 1; i <= n; i++) {
            prod *= x -i;
        }
        return logGamma1p(x - (n + 1)) + std::log(prod);
    } else {
        double sum = lanczos(x);
        double tmp = x + 607.0 / 128.0 + 0.5;
        ret = ((x + 0.5) * std::log(tmp) - tmp + HALF_LOG_2_PI + std::log(sum / x));
    }
    return ret;
}

double DigammaCache::invGamma1pm1(double x) {
    if(x < -0.5) {
        throw std::invalid_argument("X is too small.");
    }

    if(x > 1.5) {
        throw std::invalid_argument("X is too large.");
    }

    double ret;
    double t = x <= 0.5 ? x : (x - 0.5) - 0.5;
    if(t < 0.0) {
        double a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
        double b = INV_GAMMA1P_M1_B8;
        b = INV_GAMMA1P_M1_B7 + t * b;
        b = INV_GAMMA1P_M1_B6 + t * b;
        b = INV_GAMMA1P_M1_B5 + t * b;
        b = INV_GAMMA1P_M1_B4 + t * b;
        b = INV_GAMMA1P_M1_B3 + t * b;
        b = INV_GAMMA1P_M1_B2 + t * b;
        b = INV_GAMMA1P_M1_B1 + t * b;
        b = 1.0 + t * b;

        double c = INV_GAMMA1P_M1_C13 + t * (a / b);
        c = INV_GAMMA1P_M1_C12 + t * c;
        c = INV_GAMMA1P_M1_C11 + t * c;
        c = INV_GAMMA1P_M1_C10 + t * c;
        c = INV_GAMMA1P_M1_C9 + t * c;
        c = INV_GAMMA1P_M1_C8 + t * c;
        c = INV_GAMMA1P_M1_C7 + t * c;
        c = INV_GAMMA1P_M1_C6 + t * c;
        c = INV_GAMMA1P_M1_C5 + t * c;
        c = INV_GAMMA1P_M1_C4 + t * c;
        c = INV_GAMMA1P_M1_C3 + t * c;
        c = INV_GAMMA1P_M1_C2 + t * c;
        c = INV_GAMMA1P_M1_C1 + t * c;
        c = INV_GAMMA1P_M1_C + t * c;
        if (x > 0.5) {
            ret = t * c / x;
        } else {
            ret = x * ((c + 0.5) + 0.5);
        }
    } else {
        double p = INV_GAMMA1P_M1_P6;
        p = INV_GAMMA1P_M1_P5 + t * p;
        p = INV_GAMMA1P_M1_P4 + t * p;
        p = INV_GAMMA1P_M1_P3 + t * p;
        p = INV_GAMMA1P_M1_P2 + t * p;
        p = INV_GAMMA1P_M1_P1 + t * p;
        p = INV_GAMMA1P_M1_P0 + t * p;

        double q = INV_GAMMA1P_M1_Q4;
        q = INV_GAMMA1P_M1_Q3 + t * q;
        q = INV_GAMMA1P_M1_Q2 + t * q;
        q = INV_GAMMA1P_M1_Q1 + t * q;
        q = 1.0 + t * q;

        double c = INV_GAMMA1P_M1_C13 + (p / q) * t;
        c = INV_GAMMA1P_M1_C12 + t * c;
        c = INV_GAMMA1P_M1_C11 + t * c;
        c = INV_GAMMA1P_M1_C10 + t * c;
        c = INV_GAMMA1P_M1_C9 + t * c;
        c = INV_GAMMA1P_M1_C8 + t * c;
        c = INV_GAMMA1P_M1_C7 + t * c;
        c = INV_GAMMA1P_M1_C6 + t * c;
        c = INV_GAMMA1P_M1_C5 + t * c;
        c = INV_GAMMA1P_M1_C4 + t * c;
        c = INV_GAMMA1P_M1_C3 + t * c;
        c = INV_GAMMA1P_M1_C2 + t * c;
        c = INV_GAMMA1P_M1_C1 + t * c;
        c = INV_GAMMA1P_M1_C0 + t * c;

        if (x > 0.5) {
            ret = (t / x) * ((c - 0.5) - 0.5);
        } else {
            ret = x * c;
        }
    }

    return ret;
}

double DigammaCache::logGamma1p(double x) {
    if(x < -0.5) {
        throw std::invalid_argument("X is too small.");
    }

    if(x > 1.5) {
        throw std::invalid_argument("X is too large.");
    }

    return -std::log1p(invGamma1pm1(x));
}

double DigammaCache::lanczos(double x) {
    double sum = 0.0;
    for(int i = 14; i > 0; --i) {
        sum += LANCZOS[i] / (x + i);
    }
    return sum + LANCZOS[0];
}

double DigammaCache::gamma(double x) {
    if ((x == std::rint(x)) && (x <= 0.0)) {
        return nan("");
    }

    double ret;
    double absX = std::abs(x);
    if (absX <= 20.0) {
        if (x >= 1.0) {
            /*
             * From the recurrence relation
             * Gamma(x) = (x - 1) * ... * (x - n) * Gamma(x - n),
             * then
             * Gamma(t) = 1 / [1 + invGamma1pm1(t - 1)],
             * where t = x - n. This means that t must satisfy
             * -0.5 <= t - 1 <= 1.5.
             */
            double prod = 1.0;
            double t = x;
            while (t > 2.5) {
                t -= 1.0;
                prod *= t;
            }
            ret = prod / (1.0 + invGamma1pm1(t - 1.0));
        } else {
            /*
             * From the recurrence relation
             * Gamma(x) = Gamma(x + n + 1) / [x * (x + 1) * ... * (x + n)]
             * then
             * Gamma(x + n + 1) = 1 / [1 + invGamma1pm1(x + n)],
             * which requires -0.5 <= x + n <= 1.5.
             */
            double prod = x;
            double t = x;
            while (t < -0.5) {
                t += 1.0;
                prod *= t;
            }
            ret = 1.0 / (prod * (1.0 + invGamma1pm1(t)));
        }
    } else {
        double y = absX + LANCZOS_G + 0.5;
        double gammaAbs = SQRT_TWO_PI / x *
                                std::pow(y, absX + 0.5) *
                std::exp(-y) * lanczos(absX);
        if (x > 0.0) {
            ret = gammaAbs;
        } else {
            /*
             * From the reflection formula
             * Gamma(x) * Gamma(1 - x) * sin(pi * x) = pi,
             * and the recurrence relation
             * Gamma(1 - x) = -x * Gamma(-x),
             * it is found
             * Gamma(x) = -pi / [x * sin(pi * x) * Gamma(-x)].
             */
            ret = -PI /
                  (x * std::sin(PI * x) * gammaAbs);
        }
    }
    return ret;
}
