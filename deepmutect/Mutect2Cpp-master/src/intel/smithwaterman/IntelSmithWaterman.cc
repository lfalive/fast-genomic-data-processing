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
#include <vector>
#include <math.h>
#include "intel/common/avx.h"
#include "IntelSmithWaterman.h"
#include "smithwaterman_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include <assert.h>
#include "intel/common/debug.h"
#include "avx2_impl.h"
#ifndef __APPLE__
  #include "avx512_impl.h"
#endif



int32_t (*g_runSWOnePairBT)(int32_t match, int32_t mismatch, int32_t open, int32_t extend,uint8_t *seq1, uint8_t *seq2, int16_t len1, int16_t len2, int8_t overhangStrategy, char *cigarArray, int32_t cigarLen, uint32_t *cigarCount, int32_t *offset);

 void  smithwaterman_initial()
{

if(is_avx512_supported())
      {
    #ifndef __APPLE__
        DBG("Using CPU-supported AVX-512 instructions");
        g_runSWOnePairBT = runSWOnePairBT_fp_avx512;

    #else
        assert(false);
    #endif
      }
      else
      {
          g_runSWOnePairBT = runSWOnePairBT_fp_avx2;
      }
      return;
}
/*
 * Class:     com_intel_gkl_smithwaterman_IntelSmithWaterman
 * Method:    alignNative
 */
int SmithWaterman_align(uint8_t * ref, int refLength, uint8_t * alt, int altLength, uint8_t *cigar, int cigarLength, int match, int mismatch, int open, int extend, uint8_t strategy)
{

    uint32_t count = 0;
    int32_t offset = 0;

    // call the low level routine
    // Sequence length should fit in 16 bytes. This is validated earlier at the Java layer.
     int32_t result = g_runSWOnePairBT(match, mismatch, open, extend,ref,  alt, (int16_t)refLength, (int16_t)altLength, strategy, (char *) cigar, cigarLength, &count, &offset);
    return offset;
}

///*
// * Class:     com_intel_gkl_smithwaterman_IntelSmithWaterman
// * Method:    doneNative
// * Signature: ()V
// */
//JNIEXPORT void JNICALL Java_com_intel_gkl_smithwaterman_IntelSmithWaterman_doneNative
//  (JNIEnv *, jclass)
//{
//
//}
