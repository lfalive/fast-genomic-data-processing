//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_ARITHMETICUTILS_H
#define MUTECT2CPP_MASTER_ARITHMETICUTILS_H


class ArithmeticUtils {
public:
    static int gcd(int p, int q);

    static long mulAndCheck(long a, long b);

private:
    static int gcdPositive(int a, int b);

    static int numberOfTrailingZeros(unsigned int a);
};


#endif //MUTECT2CPP_MASTER_ARITHMETICUTILS_H
