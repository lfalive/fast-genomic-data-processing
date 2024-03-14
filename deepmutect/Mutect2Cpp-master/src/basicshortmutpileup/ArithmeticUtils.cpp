//
// Created by cluster on 22-11-8.
//

#include "ArithmeticUtils.h"
#include <climits>
#include <stdexcept>
#include <cmath>


int ArithmeticUtils::gcd(int p, int q) {
    int a = p;
    int b = q;
    if (a == 0 ||
        b == 0) {
        if (a == INT_MIN ||
            b == INT_MIN) {
            throw std::invalid_argument("");
        }
        return std::abs(a + b);
    }

    long al = a;
    long bl = b;
    bool useLong = false;
    if (a < 0) {
        if(INT_MIN == a) {
            useLong = true;
        } else {
            a = -a;
        }
        al = -al;
    }
    if (b < 0) {
        if (INT_MIN == b) {
            useLong = true;
        } else {
            b = -b;
        }
        bl = -bl;
    }
    if (useLong) {
        if(al == bl) {
            throw std::invalid_argument("");
        }
        long blbu = bl;
        bl = al;
        al = blbu % al;
        if (al == 0) {
            if (bl > INT_MAX) {
                throw std::invalid_argument("");
            }
            return (int) bl;
        }
        blbu = bl;

        // Now "al" and "bl" fit in an "int".
        b = (int) al;
        a = (int) (blbu % al);
    }

    return gcdPositive(a, b);
}

int ArithmeticUtils::gcdPositive(int a, int b) {
    if (a == 0) {
        return b;
    }
    else if (b == 0) {
        return a;
    }

    // Make "a" and "b" odd, keeping track of common power of 2.
    int aTwos = numberOfTrailingZeros(a);
    a >>= aTwos;
    int bTwos = numberOfTrailingZeros(b);
    b >>= bTwos;
    int shift = std::min(aTwos, bTwos);

    // "a" and "b" are positive.
    // If a > b then "gdc(a, b)" is equal to "gcd(a - b, b)".
    // If a < b then "gcd(a, b)" is equal to "gcd(b - a, a)".
    // Hence, in the successive iterations:
    //  "a" becomes the absolute difference of the current values,
    //  "b" becomes the minimum of the current values.
    while (a != b) {
        int delta = a - b;
        b = std::min(a, b);
        a = std::abs(delta);

        // Remove any power of 2 in "a" ("b" is guaranteed to be odd).
        a >>= numberOfTrailingZeros(a);
    }

    // Recover the common power of 2.
    return a << shift;
}

 int ArithmeticUtils::numberOfTrailingZeros(unsigned int i) {
    i = ~i & (i - 1);
    if (i <= 0) return i & 32;
    int n = 1;
    if (i > 1 << 16) { n += 16; i >>= 16; }
    if (i > 1 <<  8) { n +=  8; i >>=  8; }
    if (i > 1 <<  4) { n +=  4; i >>=  4; }
    if (i > 1 <<  2) { n +=  2; i >>=  2; }
    return n + ((int)i >> 1);
}

long ArithmeticUtils::mulAndCheck(long a, long b) {
    long ret;
    if (a > b) {
        // use symmetry to reduce boundary cases
        ret = mulAndCheck(b, a);
    } else {
        if (a < 0) {
            if (b < 0) {
                // check for positive overflow with negative a, negative b
                if (a >= LONG_MAX / b) {
                    ret = a * b;
                } else {
                    throw std::invalid_argument("");
                }
            } else if (b > 0) {
                // check for negative overflow with negative a, positive b
                if (LONG_MIN / b <= a) {
                    ret = a * b;
                } else {
                    throw std::invalid_argument("");

                }
            } else {
                // assert b == 0
                ret = 0;
            }
        } else if (a > 0) {
            // assert a > 0
            // assert b > 0

            // check for positive overflow with positive a, positive b
            if (a <= LONG_MAX / b) {
                ret = a * b;
            } else {
                throw std::invalid_argument("");
            }
        } else {
            // assert a == 0
            ret = 0;
        }
    }
    return ret;
}
