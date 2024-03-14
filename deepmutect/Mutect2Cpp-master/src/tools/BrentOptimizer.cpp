//
// Created by cluster on 22-11-9.
//

#include "BrentOptimizer.h"
#include <cmath>

const double BrentOptimizer::GOLDEN_SECTION = 0.5 * (3 - std::sqrt(5));

BrentOptimizer::BrentOptimizer(double rel, double abs, double min, double max, double guess, double (*fun)(double, double, const std::vector<EStep> &), int maxEval, double artifactBetaMean, const std::vector<EStep> & potentialArtifacts) : relativeThreshold(rel),
absoluteThreshold(abs), min(min), max(max), start(guess),fun(fun), maxEval(maxEval), artifactBetaMean(artifactBetaMean), potentialArtifacts(potentialArtifacts){

}

UnivariatePointValuePair BrentOptimizer::doOptimize() {
    bool isMinim = true;
    double lo = min;
    double mid = start;
    double hi = max;

    double a;
    double b;
    if (lo < hi) {
        a = lo;
        b = hi;
    } else {
        a = hi;
        b = lo;
    }

    double x = mid;
    double v = x;
    double w = x;
    double d = 0;
    double e = 0;
    double fx = computeObjectiveValue(x);
    if (!isMinim) {
        fx = -fx;
    }
    double fv = fx;
    double fw = fx;

    UnivariatePointValuePair previous = UnivariatePointValuePair::default_pair;
    UnivariatePointValuePair current
            = UnivariatePointValuePair(x, isMinim ? fx : -fx);
    // Best point encountered so far (which is the initial guess).
    UnivariatePointValuePair bestone = current;

    while (true) {
        double m = 0.5 * (a + b);
        double tol1 = relativeThreshold * std::abs(x) + absoluteThreshold;
        double tol2 = 2 * tol1;

        // Default stopping criterion.
        bool stop = std::abs(x - m) <= tol2 - 0.5 * (b - a);
        if (!stop) {
            double p = 0;
            double q = 0;
            double r = 0;
            double u = 0;

            if (std::abs(e) > tol1) { // Fit parabola.
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2 * (q - r);

                if (q > 0) {
                    p = -p;
                } else {
                    q = -q;
                }

                r = e;
                e = d;

                if (p > q * (a - x) &&
                    p < q * (b - x) &&
                    std::abs(p) < std::abs(0.5 * q * r)) {
                    // Parabolic interpolation step.
                    d = p / q;
                    u = x + d;

                    // f must not be evaluated too close to a or b.
                    if (u - a < tol2 || b - u < tol2) {
                        if (x <= m) {
                            d = tol1;
                        } else {
                            d = -tol1;
                        }
                    }
                } else {
                    // Golden section step.
                    if (x < m) {
                        e = b - x;
                    } else {
                        e = a - x;
                    }
                    d = GOLDEN_SECTION * e;
                }
            } else {
                // Golden section step.
                if (x < m) {
                    e = b - x;
                } else {
                    e = a - x;
                }
                d = GOLDEN_SECTION * e;
            }

            // Update by at least "tol1".
            if (std::abs(d) < tol1) {
                if (d >= 0) {
                    u = x + tol1;
                } else {
                    u = x - tol1;
                }
            } else {
                u = x + d;
            }

            double fu = computeObjectiveValue(u);
            if (!isMinim) {
                fu = -fu;
            }

            // User-defined convergence checker.
            previous = current;
            current = UnivariatePointValuePair(u, isMinim ? fu : -fu);
            UnivariatePointValuePair tmp = best(&previous,
                                                &current,
                                                isMinim);
            bestone = best(&bestone,
                        &tmp,
                        isMinim);


            // Update a, b, v, w and x.
            if (fu <= fx) {
                if (u < x) {
                    b = x;
                } else {
                    a = x;
                }
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            } else {
                if (u < x) {
                    a = u;
                } else {
                    b = u;
                }
                if (fu <= fw ||
                        AreDoubleSame(w, x)) {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                } else if (fu <= fv ||
                           AreDoubleSame(v, x) ||
                           AreDoubleSame(v, w)) {
                    v = u;
                    fv = fu;
                }
            }
        } else { // Default termination (Brent's criterion).
            UnivariatePointValuePair tmp = best(&previous,
                                                &current,
                                                isMinim);
            return best(&bestone,
                        &tmp,
                        isMinim);
        }
        if(currentEval > maxEval) {
            return bestone;
        }
    }
}

double BrentOptimizer::computeObjectiveValue(double x) {
    currentEval++;
    return fun(x, artifactBetaMean, potentialArtifacts);
}

UnivariatePointValuePair BrentOptimizer::best(UnivariatePointValuePair* a, UnivariatePointValuePair* b, bool isMinim) {
    if (a == nullptr) {
        return *b;
    }
    if (b == nullptr) {
        return *a;
    }

    if (isMinim) {
        return a->getValue() <= b->getValue() ? *a : *b;
    } else {
        return a->getValue() >= b->getValue() ? *a : *b;
    }
}

bool BrentOptimizer::AreDoubleSame(double dFirstVal, double dSecondVal) {
    return std::fabs(dFirstVal - dSecondVal) < std::numeric_limits<double>::epsilon();
}
