//
// Created by cluster on 22-11-8.
//

#include "CombinatoricsUtils.h"
#include "ArithmeticUtils.h"
#include <cmath>

double CombinatoricsUtils::binomialCoefficientLog(int n, int k) {
    if ((n == k) || (k == 0)) {
        return 0;
    }
    if ((k == 1) || (k == n - 1)) {
        return std::log(n);
    }

    /*
     * For values small enough to do exact integer computation,
     * return the log of the exact value
     */
    if (n < 67) {
        return std::log(binomialCoefficient(n,k));
    }

    /*
     * Return the log of binomialCoefficientDouble for values that will not
     * overflow binomialCoefficientDouble
     */
    if (n < 1030) {
        return std::log(binomialCoefficientDouble(n, k));
    }

    if (k > n / 2) {
        return binomialCoefficientLog(n, n - k);
    }

    /*
     * Sum logs for values that could overflow
     */
    double logSum = 0;

    // n!/(n-k)!
    for (int i = n - k + 1; i <= n; i++) {
        logSum += std::log(i);
    }

    // divide by k!
    for (int i = 2; i <= k; i++) {
        logSum -= std::log(i);
    }

    return logSum;
}

long CombinatoricsUtils::binomialCoefficient(int n, int k) {
    if ((n == k) || (k == 0)) {
        return 1;
    }
    if ((k == 1) || (k == n - 1)) {
        return n;
    }
    // Use symmetry for large k
    if (k > n / 2) {
        return binomialCoefficient(n, n - k);
    }

    // We use the formula
    // (n choose k) = n! / (n-k)! / k!
    // (n choose k) == ((n-k+1)*...*n) / (1*...*k)
    // which could be written
    // (n choose k) == (n-1 choose k-1) * n / k
    long result = 1;
    if (n <= 61) {
        // For n <= 61, the naive implementation cannot overflow.
        int i = n - k + 1;
        for (int j = 1; j <= k; j++) {
            result = result * i / j;
            i++;
        }
    } else if (n <= 66) {
        // For n > 61 but n <= 66, the result cannot overflow,
        // but we must take care not to overflow intermediate values.
        int i = n - k + 1;
        for (int j = 1; j <= k; j++) {
            // We know that (result * i) is divisible by j,
            // but (result * i) may overflow, so we split j:
            // Filter out the gcd, d, so j/d and i/d are integer.
            // result is divisible by (j/d) because (j/d)
            // is relative prime to (i/d) and is a divisor of
            // result * (i/d).
            long d = ArithmeticUtils::gcd(i, j);
            result = (result / (j / d)) * (i / d);
            i++;
        }
    } else {
        // For n > 66, a result overflow might occur, so we check
        // the multiplication, taking care to not overflow
        // unnecessary.
        int i = n - k + 1;
        for (int j = 1; j <= k; j++) {
            long d = ArithmeticUtils::gcd(i, j);
            result = ArithmeticUtils::mulAndCheck(result / (j / d), i / d);
            i++;
        }
    }
    return result;
}

double CombinatoricsUtils::binomialCoefficientDouble(int n, int k) {
    if ((n == k) || (k == 0)) {
        return 1;
    }
    if ((k == 1) || (k == n - 1)) {
        return n;
    }
    if (k > n/2) {
        return binomialCoefficientDouble(n, n - k);
    }
    if (n < 67) {
        return binomialCoefficient(n,k);
    }

    double result = 1;
    for (int i = 1; i <= k; i++) {
        result *= (double)(n - k + i) / (double)i;
    }

    return std::floor(result + 0.5);
}
