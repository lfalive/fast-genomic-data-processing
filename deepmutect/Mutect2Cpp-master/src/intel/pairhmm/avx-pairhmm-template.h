/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2021 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <iostream>
#include "Haplotype.h"
#include "trie/buildTreeUtils.h"

#ifdef PRECISION

void CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(const testcase& tc, int COLS, int numMaskVecs, MASK_TYPE (*maskArr)[NUM_DISTINCT_CHARS]) {

    const int maskBitCnt = MAIN_TYPE_SIZE ;

    for (int vi=0; vi < numMaskVecs; ++vi) {
        for (int rs=0; rs < NUM_DISTINCT_CHARS; ++rs) {
            maskArr[vi][rs] = 0 ;
        }
        maskArr[vi][AMBIG_CHAR] = MASK_ALL_ONES ;
    }

    for (int col=1; col < COLS; ++col) {
        int mIndex = (col-1) / maskBitCnt ;
        int mOffset = (col-1) % maskBitCnt ;
        MASK_TYPE bitMask = ((MASK_TYPE)0x1) << (maskBitCnt-1-mOffset) ;

        char hapChar = ConvertChar::get(tc.hap[col-1]);

        if (hapChar == AMBIG_CHAR) {
            for (int ci=0; ci < NUM_DISTINCT_CHARS; ++ci)
                maskArr[mIndex][ci] |= bitMask ;
        }

        maskArr[mIndex][hapChar] |= bitMask ;
        // bit corresponding to col 1 will be the MSB of the mask 0
        // bit corresponding to col 2 will be the MSB-1 of the mask 0
        // ...
        // bit corresponding to col 32 will be the LSB of the mask 0
        // bit corresponding to col 33 will be the MSB of the mask 1
        // ...
    }

}

void CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(const testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess) {

    for (int ri=0; ri < numRowsToProcess; ++ri) {
        rsArr[ri] = ConvertChar::get(tc.readForPairHmm->rs[ri+beginRowIndex-1]) ;
    }

    for (int ei=0; ei < AVX_LENGTH; ++ei) {
        lastMaskShiftOut[ei] = 0 ;
    }
}

#define SET_MASK_WORD(__dstMask, __srcMask, __lastShiftOut, __shiftBy, __maskBitCnt){ \
    MASK_TYPE __bitMask = (((MASK_TYPE)0x1) << __shiftBy) - 1 ;            \
    MASK_TYPE __nextShiftOut = (__srcMask & __bitMask) << (__maskBitCnt - __shiftBy) ; \
    __dstMask = (__srcMask >> __shiftBy) | __lastShiftOut ;        \
    __lastShiftOut = __nextShiftOut ;                    \
}


void CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)(int maskIndex, BITMASK_VEC& bitMaskVec, MASK_TYPE (*maskArr) [NUM_DISTINCT_CHARS], char* rsArr, MASK_TYPE* lastMaskShiftOut, int maskBitCnt) {

    for (int ei=0; ei < AVX_LENGTH/2; ++ei) {
        SET_MASK_WORD(bitMaskVec.getLowEntry(ei), maskArr[maskIndex][rsArr[ei]],
                lastMaskShiftOut[ei], ei, maskBitCnt) ;

        int ei2 = ei + AVX_LENGTH/2 ; // the second entry index
        SET_MASK_WORD(bitMaskVec.getHighEntry(ei), maskArr[maskIndex][rsArr[ei2]],
                lastMaskShiftOut[ei2], ei2, maskBitCnt) ;
    }

}


inline void CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (BITMASK_VEC& bitMaskVec, SIMD_TYPE& distm, SIMD_TYPE& _1_distm, SIMD_TYPE& distmChosen) {

    VEC_BLENDV(distmChosen, distm, _1_distm, bitMaskVec.getCombinedMask());

    bitMaskVec.shift_left_1bit() ;
}

/*
 * This function:
 * 1- Intializes probability values p_MM, p_XX, P_YY, p_MX, p_GAPM and pack them into vectors
 * 2- Precompute parts of "distm" which only depeneds on a row number and pack it into vector
 */

template<class NUMBER> void CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)(int ROWS, int COLS, NUMBER* shiftOutM, NUMBER *shiftOutX, NUMBER *shiftOutY, Context<NUMBER> ctx, testcase *tc)
{
    NUMBER zero = ctx._(0.0);
    // Casting is fine because the algorithm is intended to have limited precision.
    NUMBER init_Y = ctx.INITIAL_CONSTANT / (NUMBER)(tc->haplen);
    for (int s=0;s<ROWS+COLS+AVX_LENGTH;s++)
    {
        shiftOutM[s] = zero;
        shiftOutX[s] = zero;
        shiftOutY[s] = init_Y;
    }
}

/*
 * This function handles pre-stripe computation:
 * 1- Retrieve probaility vectors from memory 
 * 2- Initialize M, X, Y vectors with all 0's (for the first stripe) and shifting the last row from previous stripe for the rest 
 */

template<class NUMBER> inline void CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(
        int stripeIdx, Context<NUMBER> ctx, testcase *tc, SIMD_TYPE &pGAPM, SIMD_TYPE &pMM, SIMD_TYPE &pMX, SIMD_TYPE &pXX, SIMD_TYPE &pMY, SIMD_TYPE &pYY,
        SIMD_TYPE &rs, UNION_TYPE &rsN, SIMD_TYPE &distm, SIMD_TYPE &_1_distm,  SIMD_TYPE *distm1D, SIMD_TYPE N_packed256, SIMD_TYPE *p_MM , SIMD_TYPE *p_GAPM ,
        SIMD_TYPE *p_MX, SIMD_TYPE *p_XX , SIMD_TYPE *p_MY, SIMD_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1,
        UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER* shiftOutX, NUMBER* shiftOutM)
{
    int i = stripeIdx;
    pGAPM = p_GAPM[i];
    pMM   = p_MM[i];
    pMX   = p_MX[i];
    pXX   = p_XX[i];
    pMY   = p_MY[i];
    pYY   = p_YY[i];

    NUMBER zero = ctx._(0.0);
    // Casting is fine because the algorithm is intended to have limited precision.
    NUMBER init_Y = ctx.INITIAL_CONSTANT / (NUMBER)(tc->haplen);
    UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
    UNION_TYPE packed3;  packed3.d = VEC_SET1_VAL(3.0);

    distm = distm1D[i];
    _1_distm = VEC_SUB(packed1.d, distm);

    distm = VEC_DIV(distm, packed3.d);

    /* initialize M_t_2, M_t_1, X_t_2, X_t_1, Y_t_2, Y_t_1 */
    M_t_2.d = VEC_SET1_VAL(zero);
    X_t_2.d = VEC_SET1_VAL(zero);

    if (i==0) {
        M_t_1.d = VEC_SET1_VAL(zero);
        X_t_1.d = VEC_SET1_VAL(zero);
        Y_t_2.d = VEC_SET_LSE(init_Y);
        Y_t_1.d = VEC_SET1_VAL(zero);
    }
    else {
        X_t_1.d = VEC_SET_LSE(shiftOutX[AVX_LENGTH]);
        M_t_1.d = VEC_SET_LSE(shiftOutM[AVX_LENGTH]);
        Y_t_2.d = VEC_SET1_VAL(zero);
        Y_t_1.d = VEC_SET1_VAL(zero);
    }
    M_t_1_y = M_t_1;
}

/*
 *  This function is the main compute kernel to compute M, X and Y
 */

inline void CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(UNION_TYPE &M_t, UNION_TYPE &X_t, UNION_TYPE &Y_t, UNION_TYPE &M_t_y,
        UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1,
        SIMD_TYPE pMM, SIMD_TYPE pGAPM, SIMD_TYPE pMX, SIMD_TYPE pXX, SIMD_TYPE pMY, SIMD_TYPE pYY, SIMD_TYPE distmSel)
{
    /* Compute M_t <= distm * (p_MM*M_t_2 + p_GAPM*X_t_2 + p_GAPM*Y_t_2) */
    M_t.d = VEC_MUL(VEC_ADD(VEC_ADD(VEC_MUL(M_t_2.d, pMM), VEC_MUL(X_t_2.d, pGAPM)), VEC_MUL(Y_t_2.d, pGAPM)), distmSel);
    //M_t.d = VEC_MUL( VEC_ADD(VEC_MUL(M_t_2.d, pMM), VEC_MUL(VEC_ADD(X_t_2.d, Y_t_2.d), pGAPM)), distmSel);

    M_t_y = M_t;

    /* Compute X_t */
    X_t.d = VEC_ADD(VEC_MUL(M_t_1.d, pMX) , VEC_MUL(X_t_1.d, pXX));

    /* Compute Y_t */
    Y_t.d = VEC_ADD(VEC_MUL(M_t_1_y.d, pMY) , VEC_MUL(Y_t_1.d, pYY));
}

/*
 * This is the main compute function. It operates on the matrix in s stripe manner.
 * The stripe height is determined by the SIMD engine type. 
 * Stripe height: "AVX float": 8, "AVX double": 4
 * For each stripe the operations are anti-diagonal based. 
 * Each anti-diagonal (M_t, Y_t, X_t) depends on the two previous anti-diagonals (M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, Y_t_1).
 * Each stripe (except the fist one) depends on the last row of the previous stripe.
 * The last stripe computation handles the addition of the last row of M and X, that's the reason for loop spliting.
 */

template<class NUMBER> NUMBER CONCAT(CONCAT(compute_full_prob_,SIMD_ENGINE), PRECISION)(testcase *tc)
{
    int ROWS = tc->readForPairHmm->rslen + 1;
    int COLS = tc->haplen + 1;
    int MAVX_COUNT = (ROWS+AVX_LENGTH-1)/AVX_LENGTH;

	/* Get initialized data */
	std::shared_ptr<NUMBER> initializedData = tc->readForPairHmm->getInitializedData<NUMBER>();

	/* Probaility arrays */
	auto* p_MM = (SIMD_TYPE*)(initializedData.get());
	auto* p_XX = (SIMD_TYPE*)(initializedData.get()) + MAVX_COUNT;
	auto* p_YY = (SIMD_TYPE*)(initializedData.get()) + 2 * MAVX_COUNT;
	auto* p_MX = (SIMD_TYPE*)(initializedData.get()) + 3 * MAVX_COUNT;
	auto* p_MY = (SIMD_TYPE*)(initializedData.get()) + 4 * MAVX_COUNT;
	auto* p_GAPM = (SIMD_TYPE*)(initializedData.get()) + 5 * MAVX_COUNT;

	/* For distm precomputation */
	auto* distm1D = (SIMD_TYPE*)(initializedData.get()) + 6 * MAVX_COUNT;

    /* Carries the values from each stripe to the next stripe */
    NUMBER shiftOutM[ROWS+COLS+AVX_LENGTH], shiftOutX[ROWS+COLS+AVX_LENGTH], shiftOutY[ROWS+COLS+AVX_LENGTH];

    /* The vector to keep the anti-diagonals of M, X, Y*/
    /* Current: M_t, X_t, Y_t */
    /* Previous: M_t_1, X_t_1, Y_t_1 */
    /* Previous to previous: M_t_2, X_t_2, Y_t_2 */ 
    UNION_TYPE  M_t, M_t_1, M_t_2, X_t, X_t_1, X_t_2, Y_t, Y_t_1, Y_t_2, M_t_y, M_t_1_y;

    /* Probality vectors */
    SIMD_TYPE pGAPM, pMM, pMX, pXX, pMY, pYY;

    struct timeval start, end;
    NUMBER result_avx2;
    Context<NUMBER> ctx;
    UNION_TYPE rs , rsN;
    HAP_TYPE hap;
    SIMD_TYPE distmSel, distmChosen ;
    SIMD_TYPE distm, _1_distm;

    int r, c;
    NUMBER zero = ctx._(0.0);
    UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
    SIMD_TYPE N_packed256 = VEC_POPCVT_CHAR('N');
    // Casting is fine because the algorithm is intended to have limited precision.
    NUMBER init_Y = ctx.INITIAL_CONSTANT / (NUMBER)(tc->haplen);
    int remainingRows = (ROWS-1) % AVX_LENGTH;
    int stripe_cnt = ((ROWS-1) / AVX_LENGTH) + (remainingRows!=0);

    const int maskBitCnt = MAIN_TYPE_SIZE ;
    const int numMaskVecs = (COLS+ROWS+maskBitCnt-1)/maskBitCnt ; // ceil function

    /* Mask precomputation for distm*/
    MASK_TYPE maskArr[numMaskVecs][NUM_DISTINCT_CHARS] ;
    CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(*tc, COLS, numMaskVecs, maskArr) ;

    char rsArr[AVX_LENGTH] ;
    MASK_TYPE lastMaskShiftOut[AVX_LENGTH] ;

    /* Precompute initialization for probabilities and shift vector*/
    CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)<NUMBER>(ROWS, COLS, shiftOutM, shiftOutX, shiftOutY,
																	 ctx, tc);

    for (int i=0;i<stripe_cnt-1;i++)
    {
        //STRIPE_INITIALIZATION
        CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM);
        CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(*tc, rsArr, lastMaskShiftOut, i*AVX_LENGTH+1, AVX_LENGTH) ;
        // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors

        BITMASK_VEC bitMaskVec ;

        for (int begin_d=1;begin_d<COLS+AVX_LENGTH;begin_d+=MAIN_TYPE_SIZE)
        {
            int numMaskBitsToProcess = std::min(MAIN_TYPE_SIZE, COLS+AVX_LENGTH-begin_d) ;
            CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)((begin_d-1)/MAIN_TYPE_SIZE, bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt) ;

            for (int mbi=0; mbi < numMaskBitsToProcess; ++mbi) {
                CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (bitMaskVec, distm, _1_distm, distmChosen) ;
                int ShiftIdx = begin_d + mbi + AVX_LENGTH;

                CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                        pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(M_t, shiftOutM[ShiftIdx], shiftOutM[begin_d+mbi]);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(X_t, shiftOutX[ShiftIdx], shiftOutX[begin_d+mbi]);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(Y_t_1, shiftOutY[ShiftIdx], shiftOutY[begin_d+mbi]);

                M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;
            }
        }
    }

    int i = stripe_cnt-1;
    {
        //STRIPE_INITIALIZATION
        CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM);

        if (remainingRows==0) remainingRows=AVX_LENGTH;
        CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(*tc, rsArr, lastMaskShiftOut, i*AVX_LENGTH+1, remainingRows) ;

        SIMD_TYPE sumM, sumX;
        sumM = VEC_SET1_VAL(zero);
        sumX = VEC_SET1_VAL(zero);

        // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors
        BITMASK_VEC bitMaskVec ;

        for (int begin_d=1;begin_d<COLS+remainingRows-1;begin_d+=MAIN_TYPE_SIZE)
        {
            int numMaskBitsToProcess = std::min(MAIN_TYPE_SIZE, COLS+remainingRows-1-begin_d) ;
            CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE),PRECISION)((begin_d-1)/MAIN_TYPE_SIZE, bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt) ;
            for (int mbi=0; mbi < numMaskBitsToProcess; ++mbi) {
                CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (bitMaskVec, distm, _1_distm, distmChosen) ;
                int ShiftIdx = begin_d + mbi +AVX_LENGTH;
                CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                        pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                sumM  = VEC_ADD(sumM, M_t.d);
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(M_t, shiftOutM[ShiftIdx]);

                sumX  = VEC_ADD(sumX, X_t.d);
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(X_t, shiftOutX[ShiftIdx]);

                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(Y_t_1, shiftOutY[ShiftIdx]);

                M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;
            }
        }
        UNION_TYPE sumMX;
        sumMX.d = VEC_ADD(sumM, sumX);
        result_avx2 = sumMX.f[remainingRows-1];
    }
    // std::cout << result_avx2 << std::endl;
    return result_avx2;
}

/*
 * trie_pairhmm
 */


void CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(const trie_testcase& tc, const std::vector<int> &COLS, const std::vector<int> &numMaskVecs, MASK_TYPE ***maskArr) {

    const int maskBitCnt = MAIN_TYPE_SIZE ;
    const int haps_size = tc.haplotypeDataArray.size();

    for(int i = 0; i < haps_size; i++) {
        for (int vi=0; vi < numMaskVecs[i]; ++vi) {
            for (int rs=0; rs < NUM_DISTINCT_CHARS; ++rs) {
                maskArr[i][vi][rs] = 0 ;
            }
            maskArr[i][vi][AMBIG_CHAR] = MASK_ALL_ONES ;
        }

        uint8_t * hap = tc.haplotypeDataArray[i].haplotypeBases;
        for (int col=1; col < COLS[i]; ++col) {
            int mIndex = (col-1) / maskBitCnt ;
            int mOffset = (col-1) % maskBitCnt ;
            MASK_TYPE bitMask = ((MASK_TYPE)0x1) << (maskBitCnt-1-mOffset) ;

            char hapChar = ConvertChar::get(hap[col-1]);

            if (hapChar == AMBIG_CHAR) {
                for (int ci=0; ci < NUM_DISTINCT_CHARS; ++ci)
                    maskArr[i][mIndex][ci] |= bitMask ;
            }

            maskArr[i][mIndex][hapChar] |= bitMask ;
            // bit corresponding to col 1 will be the MSB of the mask 0
            // bit corresponding to col 2 will be the MSB-1 of the mask 0
            // ...
            // bit corresponding to col 32 will be the LSB of the mask 0
            // bit corresponding to col 33 will be the MSB of the mask 1
            // ...
        }
    }
}

template<class NUMBER> void CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)(int ROWS, const std::vector<int> &COLS, NUMBER** shiftOutM, NUMBER** shiftOutX, NUMBER** shiftOutY, Context<NUMBER> ctx, trie_testcase *tc)
{
    for(int i = 0; i < tc->haplotypeDataArray.size(); i++) {
        NUMBER zero = ctx._(0.0);
        // Casting is fine because the algorithm is intended to have limited precision.
        NUMBER init_Y = ctx.INITIAL_CONSTANT / (NUMBER)(COLS[i]-1);
        for (int s=0;s<ROWS+COLS[i]+AVX_LENGTH;s++)
        {
            shiftOutM[i][s] = zero;
            shiftOutX[i][s] = zero;
            shiftOutY[i][s] = init_Y;
        }
    }
}

template<class NUMBER> inline void CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(
        int stripeIdx, Context<NUMBER> ctx, trie_testcase *tc, SIMD_TYPE &pGAPM, SIMD_TYPE &pMM, SIMD_TYPE &pMX, SIMD_TYPE &pXX, SIMD_TYPE &pMY, SIMD_TYPE &pYY,
        SIMD_TYPE &rs, UNION_TYPE &rsN, SIMD_TYPE &distm, SIMD_TYPE &_1_distm,  SIMD_TYPE *distm1D, SIMD_TYPE N_packed256, SIMD_TYPE *p_MM , SIMD_TYPE *p_GAPM ,
        SIMD_TYPE *p_MX, SIMD_TYPE *p_XX , SIMD_TYPE *p_MY, SIMD_TYPE *p_YY, UNION_TYPE &M_t_2, UNION_TYPE &X_t_2, UNION_TYPE &M_t_1, UNION_TYPE &X_t_1,
        UNION_TYPE &Y_t_2, UNION_TYPE &Y_t_1, UNION_TYPE &M_t_1_y, NUMBER** shiftOutX, NUMBER** shiftOutM, trieNode *node)
{
    int i = stripeIdx;
    pGAPM = p_GAPM[i];
    pMM   = p_MM[i];
    pMX   = p_MX[i];
    pXX   = p_XX[i];
    pMY   = p_MY[i];
    pYY   = p_YY[i];

    NUMBER zero = ctx._(0.0);
    // Casting is fine because the algorithm is intended to have limited precision.
    int tmp = tc->haplotypeDataArray[node->getIndex()[0]].length;
    NUMBER init_Y = ctx.INITIAL_CONSTANT / ((NUMBER)tmp);
    UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
    UNION_TYPE packed3;  packed3.d = VEC_SET1_VAL(3.0);

    distm = distm1D[i];
    _1_distm = VEC_SUB(packed1.d, distm);

    distm = VEC_DIV(distm, packed3.d);

    /* initialize M_t_2, M_t_1, X_t_2, X_t_1, Y_t_2, Y_t_1 */
    M_t_2.d = VEC_SET1_VAL(zero);
    X_t_2.d = VEC_SET1_VAL(zero);
    tmp = node->getIndex()[0];

    if (i==0) {
        M_t_1.d = VEC_SET1_VAL(zero);
        X_t_1.d = VEC_SET1_VAL(zero);
        Y_t_2.d = VEC_SET_LSE(init_Y);
        Y_t_1.d = VEC_SET1_VAL(zero);
    }
    else {
        X_t_1.d = VEC_SET_LSE(shiftOutX[tmp][AVX_LENGTH]);
        M_t_1.d = VEC_SET_LSE(shiftOutM[tmp][AVX_LENGTH]);
        Y_t_2.d = VEC_SET1_VAL(zero);
        Y_t_1.d = VEC_SET1_VAL(zero);
    }
    M_t_1_y = M_t_1;
}

void CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(const trie_testcase& tc, char* rsArr, MASK_TYPE* lastMaskShiftOut, int beginRowIndex, int numRowsToProcess) {

    for (int ri=0; ri < numRowsToProcess; ++ri) {
        rsArr[ri] = ConvertChar::get(tc.readForPairHmm->rs[ri+beginRowIndex-1]) ;
    }

    for (int ei=0; ei < AVX_LENGTH; ++ei) {
        lastMaskShiftOut[ei] = 0 ;
    }
}

void CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)(int maskIndex, BITMASK_VEC& bitMaskVec, MASK_TYPE*** maskArr, char* rsArr, MASK_TYPE* lastMaskShiftOut, int maskBitCnt, trieNode *node) {
    int j = node->getIndex()[0];
    for (int ei=0; ei < AVX_LENGTH/2; ++ei) {
        SET_MASK_WORD(bitMaskVec.getLowEntry(ei), maskArr[j][maskIndex][rsArr[ei]],
                      lastMaskShiftOut[ei], ei, maskBitCnt) ;

        int ei2 = ei + AVX_LENGTH/2 ; // the second entry index
        SET_MASK_WORD(bitMaskVec.getHighEntry(ei), maskArr[j][maskIndex][rsArr[ei2]],
                      lastMaskShiftOut[ei2], ei2, maskBitCnt) ;
    }

}


template<class NUMBER> void CONCAT(CONCAT(compute_full_prob_with_trie_,SIMD_ENGINE), PRECISION)(trieNode *node, bool isLastStrip, int begin_d, std::vector<int> & COLS,
        BITMASK_VEC bitMaskVec, MASK_TYPE*** maskArr, char* rsArr, MASK_TYPE* lastMaskShiftOut, int maskBitCnt, SIMD_TYPE distm, SIMD_TYPE _1_distm, SIMD_TYPE distmChosen, UNION_TYPE M_t, UNION_TYPE X_t, UNION_TYPE Y_t, UNION_TYPE M_t_y,
                                                                                                   UNION_TYPE M_t_2, UNION_TYPE X_t_2, UNION_TYPE Y_t_2, UNION_TYPE M_t_1, UNION_TYPE X_t_1, UNION_TYPE M_t_1_y, UNION_TYPE Y_t_1,
                                                                                                   SIMD_TYPE pMM, SIMD_TYPE pGAPM, SIMD_TYPE pMX, SIMD_TYPE pXX, SIMD_TYPE pMY, SIMD_TYPE pYY, NUMBER** shiftOutX, NUMBER** shiftOutM,
                                                                                                   NUMBER** shiftOutY, SIMD_TYPE sumM, SIMD_TYPE sumX, int remainingRows, std::vector<NUMBER> & result) {
    if(!isLastStrip) {
        while(true) {
            std::vector<int>& indexs = node->getIndex();
            int numMaskBitsToProcess = std::min(MAIN_TYPE_SIZE, COLS[indexs[0]]+AVX_LENGTH-begin_d) ;
            CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE), PRECISION)((begin_d-1)/MAIN_TYPE_SIZE, bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt, node) ;
            int j = indexs[0];
            for (int mbi=0; mbi < numMaskBitsToProcess; ++mbi) {
                CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (bitMaskVec, distm, _1_distm, distmChosen) ;
                int ShiftIdx = begin_d + mbi + AVX_LENGTH;
                CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                                                                  pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(M_t, shiftOutM[j][ShiftIdx], shiftOutM[j][begin_d+mbi]);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(X_t, shiftOutX[j][ShiftIdx], shiftOutX[j][begin_d+mbi]);

                CONCAT(CONCAT(_vector_shift,SIMD_ENGINE), PRECISION)(Y_t_1, shiftOutY[j][ShiftIdx], shiftOutY[j][begin_d+mbi]);

                M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;
            }
            begin_d += MAIN_TYPE_SIZE;
            if(node->getChild().size() > 1) {
                int size = node->getChild().size();
                std::vector<trieNode *> nodes = node->getChild();
                for(int k = 1; k < size; k++) {
                    MASK_TYPE lastMaskShiftOut_tmp[AVX_LENGTH] ;
                    memcpy(lastMaskShiftOut_tmp, lastMaskShiftOut, AVX_LENGTH * sizeof(MASK_TYPE));
//                    for(int i = 0; i < AVX_LENGTH; i++) {
//                        lastMaskShiftOut_tmp[i] = lastMaskShiftOut[i];
//                    }
                    CONCAT(CONCAT(compute_full_prob_with_trie_,SIMD_ENGINE), PRECISION)(nodes[k], isLastStrip, begin_d, COLS,  bitMaskVec, maskArr, rsArr, lastMaskShiftOut_tmp, maskBitCnt, distm, _1_distm,
                                                                                           distmChosen, M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                                                                                           pMM, pGAPM, pMX, pXX, pMY, pYY, shiftOutX, shiftOutM, shiftOutY, sumM, sumX, remainingRows, result);
                }
            }
            if(!node->getChild().empty()) {
                node = node->getChild()[0];
            } else {
                if(COLS[indexs[0]]+AVX_LENGTH-begin_d < 0) {
                    break;
                }
            }
        }
    } else {
        while(true) {
            std::vector<int>& indexs = node->getIndex();
            int numMaskBitsToProcess = std::min(MAIN_TYPE_SIZE, COLS[indexs[0]]+remainingRows-1-begin_d) ;
            CONCAT(CONCAT(update_masks_for_cols_,SIMD_ENGINE),PRECISION)((begin_d-1)/MAIN_TYPE_SIZE, bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt, node) ;
            for (int mbi=0; mbi < numMaskBitsToProcess; ++mbi) {
                CONCAT(CONCAT(computeDistVec,SIMD_ENGINE), PRECISION) (bitMaskVec, distm, _1_distm, distmChosen) ;
                int ShiftIdx = begin_d + mbi +AVX_LENGTH;
                CONCAT(CONCAT(computeMXY,SIMD_ENGINE), PRECISION)(M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                                                                  pMM, pGAPM, pMX, pXX, pMY, pYY, distmChosen);

                sumM  = VEC_ADD(sumM, M_t.d);
                int j = indexs[0];
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(M_t, shiftOutM[j][ShiftIdx]);

                sumX  = VEC_ADD(sumX, X_t.d);
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(X_t, shiftOutX[j][ShiftIdx]);
                CONCAT(CONCAT(_vector_shift_last,SIMD_ENGINE), PRECISION)(Y_t_1, shiftOutY[j][ShiftIdx]);
                M_t_2 = M_t_1; M_t_1 = M_t; X_t_2 = X_t_1; X_t_1 = X_t;
                Y_t_2 = Y_t_1; Y_t_1 = Y_t; M_t_1_y = M_t_y;
            }
            begin_d += MAIN_TYPE_SIZE;
            if(node->getChild().size() > 1) {
                int size = node->getChild().size();
                std::vector<trieNode *> nodes = node->getChild();
                for(int k = 1; k < size; k++) {
                    MASK_TYPE lastMaskShiftOut_tmp[AVX_LENGTH] ;
                    memcpy(lastMaskShiftOut_tmp, lastMaskShiftOut, AVX_LENGTH * sizeof(MASK_TYPE));
//                    for(int i = 0; i < AVX_LENGTH; i++) {
//                        lastMaskShiftOut_tmp[i] = lastMaskShiftOut[i];
//                    }
                    CONCAT(CONCAT(compute_full_prob_with_trie_,SIMD_ENGINE), PRECISION)(nodes[k], isLastStrip, begin_d, COLS,  bitMaskVec, maskArr, rsArr, lastMaskShiftOut_tmp, maskBitCnt, distm, _1_distm,
                                                                                           distmChosen, M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                                                                                           pMM, pGAPM, pMX, pXX, pMY, pYY, shiftOutX, shiftOutM, shiftOutY, sumM, sumX, remainingRows, result);
                }
            }
            if(!node->getChild().empty()) {
                node = node->getChild()[0];
            } else {
                if(COLS[indexs[0]]+remainingRows-1-begin_d < 0) {
                    break;
                }
            }
        }
        UNION_TYPE sumMX;
        NUMBER result_avx2;
        sumMX.d = VEC_ADD(sumM, sumX);
        result_avx2 = sumMX.f[remainingRows-1];
        // std::cout << "node : " << node->getIndex()[0] << "  val : " << result_avx2 << std::endl;
        result[node->getIndex()[0]] = result_avx2;
    }
}

template<class NUMBER> std::vector<NUMBER> CONCAT(CONCAT(compute_full_prob_t_,SIMD_ENGINE), PRECISION)(trie_testcase *tc) {
    int ROWS = tc->readForPairHmm->rslen + 1;
    std::vector<int> COLS;
    int haps_num = tc->haplotypeDataArray.size();
    std::vector<NUMBER> result = std::vector<NUMBER>(haps_num);
	trieNode *root = tc->root;

    for(int i = 0; i < haps_num; i++) {
        // int tmp = haps[i]->getBasesLength();
        COLS.template emplace_back(tc->haplotypeDataArray[i].length + 1);
    }
    int MAVX_COUNT = (ROWS+AVX_LENGTH-1)/AVX_LENGTH;

    /* Get initialized data */
    NUMBER* initializedData = tc->readForPairHmm->getInitializedData<NUMBER>().get();

    /* Probaility arrays */
    auto* p_MM = (SIMD_TYPE*)initializedData;
    auto* p_XX = (SIMD_TYPE*)initializedData + MAVX_COUNT;
    auto* p_YY = (SIMD_TYPE*)initializedData + 2 * MAVX_COUNT;
    auto* p_MX = (SIMD_TYPE*)initializedData + 3 * MAVX_COUNT;
    auto* p_MY = (SIMD_TYPE*)initializedData + 4 * MAVX_COUNT;
    auto* p_GAPM = (SIMD_TYPE*)initializedData + 5 * MAVX_COUNT;

    /* For distm precomputation */
    auto* distm1D = (SIMD_TYPE*)initializedData + 6 * MAVX_COUNT;

    /* Carries the values from each stripe to the next stripe */
	NUMBER* shiftOutM[haps_num];
	NUMBER* shiftOutX[haps_num];
	NUMBER* shiftOutY[haps_num];

	for(int i = 0; i < haps_num; i++) {
		shiftOutM[i] = new NUMBER[ROWS+COLS[i]+AVX_LENGTH];
		shiftOutX[i] = new NUMBER[ROWS+COLS[i]+AVX_LENGTH];
		shiftOutY[i] = new NUMBER[ROWS+COLS[i]+AVX_LENGTH];
	}

    /* The vector to keep the anti-diagonals of M, X, Y*/
    /* Current: M_t, X_t, Y_t */
    /* Previous: M_t_1, X_t_1, Y_t_1 */
    /* Previous to previous: M_t_2, X_t_2, Y_t_2 */
    UNION_TYPE  M_t, M_t_1, M_t_2, X_t, X_t_1, X_t_2, Y_t, Y_t_1, Y_t_2, M_t_y, M_t_1_y;

    /* Probality vectors */
    SIMD_TYPE pGAPM, pMM, pMX, pXX, pMY, pYY;

    struct timeval start, end;
    Context<NUMBER> ctx;
    UNION_TYPE rs , rsN;
    HAP_TYPE hap;
    SIMD_TYPE distmSel, distmChosen ;
    SIMD_TYPE distm, _1_distm;

    int r, c;
    NUMBER zero = ctx._(0.0);
    UNION_TYPE packed1;  packed1.d = VEC_SET1_VAL(1.0);
    SIMD_TYPE N_packed256 = VEC_POPCVT_CHAR('N');
    // Casting is fine because the algorithm is intended to have limited precision.
    std::vector<NUMBER> init_Y;
    for(int i = 0; i < haps_num; i++) {
        init_Y.template emplace_back(ctx.INITIAL_CONSTANT / (NUMBER)(COLS[i]-1));
    }
    int remainingRows = (ROWS-1) % AVX_LENGTH;
    int stripe_cnt = ((ROWS-1) / AVX_LENGTH) + (remainingRows!=0);

    const int maskBitCnt = MAIN_TYPE_SIZE ;
    std::vector<int> numMaskVecs;
    for(int i = 0; i < haps_num; i++) {
        numMaskVecs.template emplace_back((COLS[i]+ROWS+maskBitCnt-1)/maskBitCnt);
    }

    /* Mask precomputation for distm*/
    MASK_TYPE *** maskArr = new MASK_TYPE**[haps_num];
    for(int i = 0; i < haps_num; i++) {
        maskArr[i] = new MASK_TYPE*[numMaskVecs[i]];
        for(int j = 0; j < numMaskVecs[i]; j++) {
            maskArr[i][j] = new MASK_TYPE[NUM_DISTINCT_CHARS];
        }
    }
    CONCAT(CONCAT(precompute_masks_,SIMD_ENGINE), PRECISION)(*tc, COLS, numMaskVecs, maskArr) ;
    char rsArr[AVX_LENGTH] ;
    MASK_TYPE lastMaskShiftOut[AVX_LENGTH] ;

    /* Precompute initialization for probabilities and shift vector*/
    CONCAT(CONCAT(initializeVectors,SIMD_ENGINE), PRECISION)<NUMBER>(ROWS, COLS, shiftOutM, shiftOutX, shiftOutY,
                                                                     ctx, tc);

    for (int i=0;i<stripe_cnt-1;i++)
    {
        for(trieNode *node : root->getChild()) {
            for(trieNode *tmp : node->getChild()) {
                CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                                                                            p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM, tmp);
                CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(*tc, rsArr, lastMaskShiftOut, i*AVX_LENGTH+1, AVX_LENGTH) ;
                // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors

                BITMASK_VEC bitMaskVec ;
                SIMD_TYPE sumM, sumX;
                sumM = VEC_SET1_VAL(zero);
                sumX = VEC_SET1_VAL(zero);
                CONCAT(CONCAT(compute_full_prob_with_trie_,SIMD_ENGINE), PRECISION)(tmp, false, 1, COLS,  bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt, distm, _1_distm,
                                                                                       distmChosen, M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                                                                                       pMM, pGAPM, pMX, pXX, pMY, pYY, shiftOutX, shiftOutM, shiftOutY, sumM, sumX, AVX_LENGTH, result);
            }
        }
    }

    int i = stripe_cnt-1;
    {
        for(trieNode *node : root->getChild()) {
            for(trieNode *tmp : node->getChild()) {
                CONCAT(CONCAT(stripeINITIALIZATION,SIMD_ENGINE), PRECISION)(i, ctx, tc, pGAPM, pMM, pMX, pXX, pMY, pYY, rs.d, rsN, distm, _1_distm, distm1D, N_packed256, p_MM , p_GAPM ,
                                                                            p_MX, p_XX , p_MY, p_YY, M_t_2, X_t_2, M_t_1, X_t_1, Y_t_2, Y_t_1, M_t_1_y, shiftOutX, shiftOutM, tmp);

                if (remainingRows==0) remainingRows=AVX_LENGTH;
                CONCAT(CONCAT(init_masks_for_row_,SIMD_ENGINE), PRECISION)(*tc, rsArr, lastMaskShiftOut, i*AVX_LENGTH+1, remainingRows) ;

                SIMD_TYPE sumM, sumX;
                sumM = VEC_SET1_VAL(zero);
                sumX = VEC_SET1_VAL(zero);

                // Since there are no shift intrinsics in AVX, keep the masks in 2 SSE vectors
                BITMASK_VEC bitMaskVec ;

                CONCAT(CONCAT(compute_full_prob_with_trie_,SIMD_ENGINE), PRECISION)(tmp, true, 1, COLS,  bitMaskVec, maskArr, rsArr, lastMaskShiftOut, maskBitCnt, distm, _1_distm,
                                                                                       distmChosen, M_t, X_t, Y_t, M_t_y, M_t_2, X_t_2, Y_t_2, M_t_1, X_t_1, M_t_1_y, Y_t_1,
                                                                                       pMM, pGAPM, pMX, pXX, pMY, pYY, shiftOutX, shiftOutM, shiftOutY, sumM, sumX, remainingRows, result);
            }
        }
    }
	for(int j = 0; j < haps_num; j++) {
		delete[] shiftOutM[j];
		delete[] shiftOutX[j];
		delete[] shiftOutY[j];
	}
	for(int k = 0; k < haps_num; k++) {
		for(int j = 0; j < numMaskVecs[k]; j++) {
			delete[] maskArr[k][j];
		}
		delete[] maskArr[k];
	}
	delete[] maskArr;
    return result;
}

#endif

