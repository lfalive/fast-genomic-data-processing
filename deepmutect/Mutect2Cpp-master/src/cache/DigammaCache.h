//
// Created by 梦想家xixi on 2021/11/1.
//

#ifndef MUTECT2CPP_MASTER_DIGAMMACACHE_H
#define MUTECT2CPP_MASTER_DIGAMMACACHE_H

#include "IntToDoubleFunctionCache.h"
#include <cmath>

class DigammaCache : public IntToDoubleFunctionCache {
private:
    static const int CACHE_SIZE = 100000;

    constexpr static const double INV_GAMMA1P_M1_A0 = .611609510448141581788E-08;

    /** The constant {@code A1} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_A1 = .624730830116465516210E-08;

    /** The constant {@code B1} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B1 = .203610414066806987300E+00;

    /** The constant {@code B2} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B2 = .266205348428949217746E-01;

    /** The constant {@code B3} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B3 = .493944979382446875238E-03;

    /** The constant {@code B4} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B4 = -.851419432440314906588E-05;

    /** The constant {@code B5} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B5 = -.643045481779353022248E-05;

    /** The constant {@code B6} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B6 = .992641840672773722196E-06;

    /** The constant {@code B7} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B7 = -.607761895722825260739E-07;

    /** The constant {@code B8} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_B8 = .195755836614639731882E-09;

    /** The constant {@code P0} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P0 = .6116095104481415817861E-08;

    /** The constant {@code P1} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P1 = .6871674113067198736152E-08;

    /** The constant {@code P2} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P2 = .6820161668496170657918E-09;

    /** The constant {@code P3} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P3 = .4686843322948848031080E-10;

    /** The constant {@code P4} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P4 = .1572833027710446286995E-11;

    /** The constant {@code P5} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P5 = -.1249441572276366213222E-12;

    /** The constant {@code P6} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_P6 = .4343529937408594255178E-14;

    /** The constant {@code Q1} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_Q1 = .3056961078365221025009E+00;

    /** The constant {@code Q2} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_Q2 = .5464213086042296536016E-01;

    /** The constant {@code Q3} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_Q3 = .4956830093825887312020E-02;

    /** The constant {@code Q4} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_Q4 = .2692369466186361192876E-03;

    /** The constant {@code C} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C = -.422784335098467139393487909917598E+00;

    /** The constant {@code C0} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C0 = .577215664901532860606512090082402E+00;

    /** The constant {@code C1} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C1 = -.655878071520253881077019515145390E+00;

    /** The constant {@code C2} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C2 = -.420026350340952355290039348754298E-01;

    /** The constant {@code C3} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C3 = .166538611382291489501700795102105E+00;

    /** The constant {@code C4} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C4 = -.421977345555443367482083012891874E-01;

    /** The constant {@code C5} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C5 = -.962197152787697356211492167234820E-02;

    /** The constant {@code C6} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C6 = .721894324666309954239501034044657E-02;

    /** The constant {@code C7} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C7 = -.116516759185906511211397108401839E-02;

    /** The constant {@code C8} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C8 = -.215241674114950972815729963053648E-03;

    /** The constant {@code C9} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C9 = .128050282388116186153198626328164E-03;

    /** The constant {@code C10} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C10 = -.201348547807882386556893914210218E-04;

    /** The constant {@code C11} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C11 = -.125049348214267065734535947383309E-05;

    /** The constant {@code C12} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C12 = .113302723198169588237412962033074E-05;

    /** The constant {@code C13} defined in {@code DGAM1}. */
    constexpr static const double INV_GAMMA1P_M1_C13 = -.205633841697760710345015413002057E-06;

    static double LANCZOS[15];


protected:
    int maxSize() override {return CACHE_SIZE;}

    double compute(int n) override;

public:
    static double digamma(double x);

    static double logGamma(double x);

    static double invGamma1pm1(double x);

    static double logGamma1p(double x);

    static double lanczos(double x);

    static double gamma(double x);

    constexpr static double LANCZOS_G = 607.0 / 128.0;

    constexpr static double SQRT_TWO_PI = 2.506628274631000502;

    constexpr static double PI = 105414357.0 / 33554432.0 + 1.984187159361080883e-9;
};


#endif //MUTECT2CPP_MASTER_DIGAMMACACHE_H
