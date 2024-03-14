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
#ifdef linux
#include <omp.h>
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <debug.h>
#include <avx.h>
#include <cassert>
#include "IntelPairHmm.h"
#include "pairhmm_common.h"
#include "avx_impl.h"
#ifndef __APPLE__
  #include "avx512_impl.h"
#endif
#include "Context.h"
#include <thread>
#include "boost/utility.hpp"
#include "utils/pairhmm/PairHMMConcurrentControl.h"
#include "trie/buildTreeUtils.h"

bool g_use_double;
int g_max_threads;

Context<float> g_ctxf;
Context<double> g_ctxd;

float (*g_compute_full_prob_float)(testcase *tc);
std::vector<float> (*g_compute_full_prob_t_float)(trie_testcase *tc);
double (*g_compute_full_prob_double)(testcase *tc);
std::vector<double> (*g_compute_full_prob_t_double)(trie_testcase *tc);

/*
 * Class:     com_intel_gkl_pairhmm_IntelPairHmm
 * Method:    initNative
 * Signature: (ZI)V
 */
/*
JNIEXPORT void JNICALL Java_com_intel_gkl_pairhmm_IntelPairHmm_initNative
(JNIEnv* env, jclass cls, jclass readDataHolder, jclass haplotypeDataHolder,
 jboolean use_double, jint max_threads)
{
  DBG("Enter");

  JavaData javaData(env);
  try {
    javaData.init(readDataHolder, haplotypeDataHolder);
  } catch (JavaException& e) {
    env->ExceptionClear();
    env->ThrowNew(env->FindClass(e.classPath), e.message);
    return;
  }

  g_use_double = use_double;

#ifdef _OPENMP
  int avail_threads = omp_get_max_threads();
  int req_threads = max_threads;
  g_max_threads = std::min(req_threads, avail_threads);

  DBG("Available threads: %d", avail_threads);
  DBG("Requested threads: %d", req_threads);
  if (req_threads > avail_threads) {
    DBG("Using %d available threads, but %d were requested", g_max_threads, req_threads);
  }
  else {
    DBG("Using %d threads", g_max_threads);
  }
#else
  if (max_threads != 1) {
    DBG("Ignoring request for %d threads; not using OpenMP implementation", max_threads);
  }
#endif


  // enable FTZ
  if (_MM_GET_FLUSH_ZERO_MODE() != _MM_FLUSH_ZERO_ON) {
    DBG("Flush-to-zero (FTZ) is enabled when running PairHMM");
  }
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  // set function pointers
  if(is_avx512_supported())
  {
#ifndef __APPLE__
    DBG("Using CPU-supported AVX-512 instructions");
    g_compute_full_prob_float = compute_fp_avx512s;
    g_compute_full_prob_double = compute_fp_avx512d;
#else
    assert(false);
#endif
  }
  else
  {
    g_compute_full_prob_float = compute_fp_avxs;
    g_compute_full_prob_double = compute_fp_avxd;
  }

  // init convert char table
  ConvertChar::init();
  DBG("Exit");
}
*/

/*
 * Class:     com_intel_gkl_pairhmm_IntelPairHmm
 * Method:    computeLikelihoodsNative
 * Signature: (II[[B[[B[[B[[B[[B[[B[D)V
 */
/*JNIEXPORT void JNICALL Java_com_intel_gkl_pairhmm_IntelPairHmm_computeLikelihoodsNative
(JNIEnv* env, jobject obj,
 jobjectArray readDataArray, jobjectArray haplotypeDataArray, jdoubleArray likelihoodArray)
{
  DBG("Enter");

  //==================================================================
  // get Java data
  JavaData javaData(env);

  std::vector<testcase> testcases;
  double* javaResults;

  try {
    testcases = javaData.getData(readDataArray, haplotypeDataArray);
    javaResults = javaData.getOutputArray(likelihoodArray);
  } catch (JavaException& e) {
    env->ExceptionClear();
    env->ThrowNew(env->FindClass(e.classPath), e.message);
    return;
  }

  //==================================================================
  // calcutate pairHMM

    try {
        #ifdef _OPENMP
              #pragma omp parallel for schedule(dynamic, 1) num_threads(g_max_threads)
        #endif
        for (int i = 0; i < testcases.size(); i++) {

            double result_final = 0;
            float result_float = g_use_double ? 0.0f : g_compute_full_prob_float(&testcases[i]);

            if (result_float < MIN_ACCEPTED) {
              double result_double = g_compute_full_prob_double(&testcases[i]);
              result_final = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT;
            }
            else {
              result_final = (double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT);
            }

            javaResults[i] = result_final;
            DBG("result = %e", result_final);
          }
    }
    catch (JavaException& e) {
       #ifdef _OPENMP
            #pragma omp barrier
       #endif
       env->ExceptionClear();
       env->ThrowNew(env->FindClass(e.classPath), "Error in pairhmm processing.");
       return;
    }

  DBG("Exit");
}*/


/*
 * Class:     com_intel_gkl_pairhmm_IntelPairHmm
 * Method:    doneNative
 * Signature: ()V
 */
//JNIEXPORT void JNICALL Java_com_intel_gkl_pairhmm_IntelPairHmm_doneNative
//(JNIEnv* env, jobject obj)
//{
//}

void initNative(bool use_double, int max_threads)
{
    g_use_double = use_double;

#ifdef _OPENMP
    int avail_threads = omp_get_max_threads();
    int req_threads = max_threads;
    g_max_threads = std::min(req_threads, avail_threads);

    DBG("Available threads: %d", avail_threads);
    DBG("Requested threads: %d", req_threads);
    if (req_threads > avail_threads) {
        DBG("Using %d available threads, but %d were requested", g_max_threads, req_threads);
    }
    else {
        DBG("Using %d threads", g_max_threads);
    }
#else
    if (max_threads != 1) {
        DBG("Ignoring request for %d threads; not using OpenMP implementation", max_threads);
    }
#endif


    // enable FTZ
    if (_MM_GET_FLUSH_ZERO_MODE() != _MM_FLUSH_ZERO_ON) {
        DBG("Flush-to-zero (FTZ) is enabled when running PairHMM");
    }
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

    // set function pointers
    if(is_avx512_supported())
    {
#ifndef __APPLE__
        DBG("Using CPU-supported AVX-512 instructions");
        g_compute_full_prob_float = compute_fp_avx512s;
        g_compute_full_prob_double = compute_fp_avx512d;
        g_compute_full_prob_t_float = compute_fp_t_avx512s;
        g_compute_full_prob_t_double = compute_fp_t_avx512d;
#else
        assert(false);
#endif
    }
    else
    {
        g_compute_full_prob_float = compute_fp_avxs;
        g_compute_full_prob_double = compute_fp_avxd;
        g_compute_full_prob_t_float = compute_fp_t_avxs;
        g_compute_full_prob_t_double = compute_fp_t_avxd;
    }

    // init convert char table
    ConvertChar::init();
    DBG("Exit");
}

// native implementation of PairHMM algorithm
void computeLikelihoodsNative(std::vector<testcase>& testcases, std::vector<double>& likelihoodArray)
{
    DBG("Enter");

    //==================================================================
    // prepare data

    //==================================================================
    // calcutate pairHMM

    try {
        #ifdef _OPENMP
              #pragma omp parallel for schedule(dynamic, 1) num_threads(g_max_threads)
        #endif
        for (int i = 0; i < testcases.size(); i++) {

            double result_final = 0;
            float result_float = g_use_double ? 0.0f : g_compute_full_prob_float(&testcases[i]);

            if (result_float < MIN_ACCEPTED) {
                double result_double = g_compute_full_prob_double(&testcases[i]);
                result_final = log10(result_double) - Context<double>::LOG10_INITIAL_CONSTANT;
            }
            else {
                result_final = (double)(log10f(result_float) - Context<float>::LOG10_INITIAL_CONSTANT);
            }

            //javaResults[i] = result_final;
            likelihoodArray[i] = result_final;
            DBG("result = %e", result_final);
        }
    } catch (const char* message) {
        std::cout << message << std::endl;
    }
}

// native implementation of PairHMM algorithm modified by hlf
void computeLikelihoodsNative_concurrent(std::vector<testcase>& testcases, std::vector<double>& likelihoodArray) {
	for (int i = 0; i < testcases.size(); i++) {
		if (BOOST_LIKELY((!PairHMMConcurrentControl::startPairHMMConcurrentMode) || testcases.size() - i < 1024)) {
			computeLikelihoodsNative_concurrent_i(testcases, likelihoodArray, i);
		} else {
			std::shared_ptr<LikelihoodsTask> likelihoods = std::make_shared<LikelihoodsTask>(testcases, likelihoodArray, (unsigned long)i, testcases.size());
			PairHMMConcurrentControl::pairHMMMutex.lock();
			PairHMMConcurrentControl::pairHMMTaskQueue.push(likelihoods);
			//std::cout << "push " + std::to_string(PairHMMConcurrentControl::pairHMMTaskQueue.size()) + '\n';
			PairHMMConcurrentControl::pairHMMMutex.unlock();
			//std::cout << std::to_string(testcases.size()) + " testcases\tpush [" + std::to_string(i) + ", " + std::to_string(testcases.size() - 1) + "] size: " +
			//                                                                                                                                          std::to_string(testcases.size() - i ) + '\n';

			for (unsigned long ind = likelihoods->index++; ind < likelihoods->testcasesSize; ind = likelihoods->index++) {
				//std::cout << std::to_string(ind) + '\n';
				computeLikelihoodsNative_concurrent_i(likelihoods->taskTestcases, likelihoods->taskLikelihoodArray, ind);
				likelihoods->count++;
			}

			// make sure all calculations have been completed
			while (likelihoods->count != likelihoods->testcasesSize)
				std::this_thread::yield();

			// Check the results of concurrent calculations
			/*std::vector<double> newArray(testcases.size());
			for (int new_i = 0; new_i < testcases.size(); new_i++) {
				computeLikelihoodsNative_concurrent_i(&testcases, &newArray, new_i);
				assert(abs(newArray[new_i] - likelihoodArray[new_i]) < 0.0000000001);
			}*/

			//std::cout << std::to_string(testcases.size()) + " testcases done.\n";
			return;
		}
	}
}

void computeLikelihoodsNative_concurrent_i(std::vector<testcase> &testcases, std::vector<double> &likelihoodArray, unsigned long i) {
//	 if ((*testcases)[i].haplen > 64 && (*testcases)[i].rslen > 64)
//		test_compute(testcases,likelihoodArray,i);

//	std::cout << "double\t" << Context<double>::INITIAL_CONSTANT << '\t' << Context<double>::LOG10_INITIAL_CONSTANT << '\t' << log10(compute_full_prob_double(&(*testcases)[i])) - Context<double>::LOG10_INITIAL_CONSTANT << std::endl;
//	std::cout << "float\t" << Context<float>::INITIAL_CONSTANT << '\t' << Context<float>::LOG10_INITIAL_CONSTANT << '\t' << (double)(log10f(compute_full_prob_float(&(*testcases)[i]))) - Context<float>::LOG10_INITIAL_CONSTANT << std::endl << std::endl;

	double result_final;
	float result_float = (BOOST_UNLIKELY(g_use_double)) ? 0.0f : g_compute_full_prob_float(&testcases[i]);

	if (BOOST_UNLIKELY(result_float < MIN_ACCEPTED)) {
		double result_double = g_compute_full_prob_double(&testcases[i]);
		result_final = log10(result_double) - Context<double>::LOG10_INITIAL_CONSTANT;
//		PairHMMConcurrentControl::compute_double_cases++;
	}
	else {
		result_final = (double)(log10f(result_float) - Context<float>::LOG10_INITIAL_CONSTANT);
	}
	likelihoodArray[i] = result_final;
}

//void test_compute(std::vector<testcase> *testcases, std::vector<double> *likelihoodArray, unsigned long i){
//	auto startTime = std::chrono::system_clock::now(), endTime = std::chrono::system_clock::now();
//
////	startTime = std::chrono::system_clock::now();
////	compute_full_prob_Fixed64(&(*testcases)[i]);
////	endTime = std::chrono::system_clock::now();
////	std::cout << "Fixed64\t" << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << std::endl;
//
//	startTime = std::chrono::system_clock::now();
//	compute_full_prob_float(&(*testcases)[i]);
//	endTime = std::chrono::system_clock::now();
//	std::cout << "float native\t" << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << std::endl;
//
//	startTime = std::chrono::system_clock::now();
//	g_compute_full_prob_float(&(*testcases)[i]);
//	endTime = std::chrono::system_clock::now();
//	std::cout << "float intel\t" << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << std::endl;
//
//	startTime = std::chrono::system_clock::now();
//	compute_full_prob_double(&(*testcases)[i]);
//	endTime = std::chrono::system_clock::now();
//	std::cout << "double native\t" << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << std::endl;
//
//	startTime = std::chrono::system_clock::now();
//	g_compute_full_prob_double(&(*testcases)[i]);
//	endTime = std::chrono::system_clock::now();
//	std::cout << "double intel\t" << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << std::endl;
//
//	std::cout << std::endl;
//}
//
//
//
//float compute_full_prob_float(testcase *tc) {
//	/*float ph2pr[128];
//	for (int i = 0; i < 128; i++)
//		ph2pr[i] = powf(10.f, -((float)i) / 10.f);
//	*/
//
//	int ROWS = tc->rslen + 1, COLS = tc->haplen + 1;
//
//	float distm[ROWS];
//	distm[0] = -1;  // distm[0] will not be used
//
//	float M[ROWS][COLS], X[ROWS][COLS], Y[ROWS][COLS], p[ROWS][6];
//	std::fill_n(p[0], 6, 0);
//	int MM = 0, GapM = 1, MX = 2, XX = 3, MY = 4, YY = 5;   //use #define ?
//	for (int r = 1; r < ROWS; r++) {
//		float delta = Context<float>::ph2pr[tc->d[r - 1] & 127];
//		float iota = Context<float>::ph2pr[tc->i[r - 1] & 127];
//		float epsilon = Context<float>::ph2pr[tc->c[r - 1] & 127];
//		p[r][MM] = 1.0f - iota - delta; // can be calculated and stored?
//		p[r][GapM] = 1.0f - epsilon;    // can be calculated and stored?
//		p[r][MX] = iota;
//		p[r][MY] = delta;
//		p[r][XX] = p[r][YY] = epsilon;
//
//		distm[r] = Context<float>::ph2pr[tc->q[r - 1] & 127];
//	}
//
//	// float init_Y = ldexpf(1.f, 120.f) / (float)tc->haplen;
//	// in double version, use ldexp(1.0, 1020.0)
//	float init_Y = Context<float>::INITIAL_CONSTANT / (float) tc->haplen;   // 防止数值下溢
//	for (int c = 0; c < COLS; c++) {
//		M[0][c] = 0;
//		X[0][c] = 0;
//		Y[0][c] = init_Y;
//	}
//
//	for (int r = 1; r < ROWS; r++) {
//		M[r][0] = 0;
//		X[r][0] = 0;
//		Y[r][0] = 0;
//	}
//
//	for (int r = 1; r < ROWS; r++) {
//		for (int c = 1; c < COLS; c++) {
//			// char _rs = tc->rs[r - 1];
//			// char _hap = tc->hap[c - 1];
//			// float _distm = _rs == _hap ? 1.0f - distm[r] : distm[r] / 3;
//
//			M[r][c] = (tc->rs[r - 1] == tc->hap[c - 1] ? 1.0f - distm[r] : distm[r] / 3)
//			          * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
//
//			X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
//
//			Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];
//		}
//	}
//
//	float result = 0;
//	for (int c = 1; c < COLS; c++)
//		result += M[ROWS - 1][c] + X[ROWS - 1][c];
//
//	/*std::cout.precision(5);
//	std::cout.setf(std::ios::scientific);
//	for (int r = 0; r < ROWS; ++r)
//	{
//		for (int c = 0; c < COLS; ++c)
//		{
//			std::cout<<M[r][c]<<"\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << "============\n";
//
//	for (int r = 0; r < ROWS; ++r)
//	{
//		for (int c = 0; c < COLS; ++c)
//		{
//			std::cout<<X[r][c]<<"\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << "============\n";
//
//	for (int r = 0; r < ROWS; ++r)
//	{
//		for (int c = 0; c < COLS; ++c)
//		{
//			std::cout<<Y[r][c]<<"\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << "============\n";*/
//
//	return result;
//}
//
//double compute_full_prob_double(testcase *tc) {
//	int ROWS = tc->rslen + 1, COLS = tc->haplen + 1;
//
//	double distm[ROWS];
//	distm[0] = -1;  // distm[0] will not be used
//
//	double M[ROWS][COLS], X[ROWS][COLS], Y[ROWS][COLS], p[ROWS][6];
//	std::fill_n(p[0], 6, 0);
//	int MM = 0, GapM = 1, MX = 2, XX = 3, MY = 4, YY = 5;   //use #define ?
//	for (int r = 1; r < ROWS; r++) {
//		double delta = Context<double>::ph2pr[tc->d[r - 1] & 127];
//		double iota = Context<double>::ph2pr[tc->i[r - 1] & 127];
//		double epsilon = Context<double>::ph2pr[tc->c[r - 1] & 127];
//		p[r][MM] = 1.0f - iota - delta;
//		p[r][GapM] = 1.0f - epsilon;
//		p[r][MX] = iota;
//		p[r][MY] = delta;
//		p[r][XX] = p[r][YY] = epsilon;
//
//		distm[r] = Context<double>::ph2pr[tc->q[r - 1] & 127];
//	}
//
//	// in double version, use ldexp(1.0, 1020.0)
//	double init_Y = Context<double>::INITIAL_CONSTANT / (double) tc->haplen;   // 防止数值下溢
//	for (int c = 0; c < COLS; c++) {
//		M[0][c] = 0;
//		X[0][c] = 0;
//		Y[0][c] = init_Y;
//	}
//
//	for (int r = 1; r < ROWS; r++) {
//		M[r][0] = 0;
//		X[r][0] = 0;
//		Y[r][0] = 0;
//	}
//
//	for (int r = 1; r < ROWS; r++) {
//		for (int c = 1; c < COLS; c++) {
//			// char _rs = tc->rs[r - 1];
//			// char _hap = tc->hap[c - 1];
//			// float _distm = _rs == _hap ? 1.0f - distm[r] : distm[r] / 3;
//
//			M[r][c] = (tc->rs[r - 1] == tc->hap[c - 1] ? 1.0f - distm[r] : distm[r] / 3)
//			          * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
//
//			X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
//
//			Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];
//		}
//	}
//
//	double result = 0;
//	for (int c = 1; c < COLS; c++)
//		result += M[ROWS - 1][c] + X[ROWS - 1][c];
//
//	/*std::cout.precision(5);
//	std::cout.setf(std::ios::scientific);
//	for (int r = 0; r < ROWS; ++r)
//	{
//		for (int c = 0; c < COLS; ++c)
//		{
//			std::cout<<M[r][c]<<"\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << "============\n";
//
//	for (int r = 0; r < ROWS; ++r)
//	{
//		for (int c = 0; c < COLS; ++c)
//		{
//			std::cout<<X[r][c]<<"\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << "============\n";
//
//	for (int r = 0; r < ROWS; ++r)
//	{
//		for (int c = 0; c < COLS; ++c)
//		{
//			std::cout<<Y[r][c]<<"\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << "============\n";*/
//
//	return result;
//}

/*
 * author: hlf
 * get Fixed64.h from https://github.com/XMunkki/FixPointCS to use this function.
 * This function has been discarded because of its low efficiency and low accuracy.
 */

/*
double compute_full_prob_Fixed64(testcase *tc) {

	int INITIAL_CONSTANT_INT = 1 << 30;
	Fixed64::FP_LONG INITIAL_CONSTANT = Fixed64::FromInt(INITIAL_CONSTANT_INT);
	float LOG10_INITIAL_CONSTANT = log10f((float)INITIAL_CONSTANT_INT);

	Fixed64::FP_LONG ph2pr[128];
	Fixed64::FP_LONG TEN = Fixed64::FromInt(10);
	for (int i = 0; i < 128; i++) {
		ph2pr[i] = Fixed64::Pow(TEN, Fixed64::Div(Fixed64::FromInt(-i), TEN));
	}

	int ROWS = tc->rslen + 1, COLS = tc->haplen + 1;

	Fixed64::FP_LONG distm[ROWS];
	distm[0] = Fixed64::MinValue;  // distm[0] will not be used

	Fixed64::FP_LONG M[ROWS][COLS], X[ROWS][COLS], Y[ROWS][COLS], p[ROWS][6];
	std::fill_n(p[0], 6, 0);

	int MM = 0, GapM = 1, MX = 2, XX = 3, MY = 4, YY = 5;   //use #define ?
	for (int r = 1; r < ROWS; r++) {
		Fixed64::FP_LONG delta = ph2pr[tc->d[r - 1] & 127];
		Fixed64::FP_LONG iota = ph2pr[tc->i[r - 1] & 127];
		Fixed64::FP_LONG epsilon = ph2pr[tc->c[r - 1] & 127];
		p[r][MM] = Fixed64::One - iota - delta; // can be calculated and stored?
		p[r][GapM] = Fixed64::One - epsilon;    // can be calculated and stored?
		p[r][MX] = iota;
		p[r][MY] = delta;
		p[r][XX] = p[r][YY] = epsilon;

		distm[r] = ph2pr[tc->q[r - 1] & 127];
	}

	// float init_Y = ldexpf(1.f, 120.f) / (float)tc->haplen;
	// in double version, use ldexp(1.0, 1020.0)
	Fixed64::FP_LONG init_Y = Fixed64::Div(INITIAL_CONSTANT, Fixed64::FromInt(tc->haplen));   // 防止数值下溢
	for (int c = 0; c < COLS; c++) {
		M[0][c] = 0;
		X[0][c] = 0;
		Y[0][c] = init_Y;
	}

	for (int r = 1; r < ROWS; r++) {
		M[r][0] = 0;
		X[r][0] = 0;
		Y[r][0] = 0;
	}

	for (int r = 1; r < ROWS; r++)
		for (int c = 1; c < COLS; c++) {
			// char _rs = tc->rs[r - 1];
			// char _hap = tc->hap[c - 1];
			// float _distm = (_rs == _hap || _rs == 'N' || _hap == 'N') ? 1.0f - distm[r] : distm[r] / 3;

			// M[r][c] = (tc->rs[r - 1] == tc->hap[c - 1]  ? 1.0f - distm[r] : distm[r] / 3)
			//           * (M[r - 1][c - 1] * p[r][MM] + X[r - 1][c - 1] * p[r][GapM] + Y[r - 1][c - 1] * p[r][GapM]);
			Fixed64::FP_LONG _distm = tc->rs[r - 1] == tc->hap[c - 1] ? Fixed64::One - distm[r] : Fixed64::Div(distm[r],Fixed64::Three);
			M[r][c] = Fixed64::Mul(_distm, Fixed64::Mul(M[r - 1][c - 1], p[r][MM]) + Fixed64::Mul(X[r - 1][c - 1] + Y[r - 1][c - 1], p[r][GapM]));

			// X[r][c] = M[r - 1][c] * p[r][MX] + X[r - 1][c] * p[r][XX];
			X[r][c] = Fixed64::Mul(M[r - 1][c], p[r][MX]) + Fixed64::Mul(X[r - 1][c], p[r][XX]);

			// Y[r][c] = M[r][c - 1] * p[r][MY] + Y[r][c - 1] * p[r][YY];
			Y[r][c] = Fixed64::Mul(M[r][c - 1], p[r][MY]) + Fixed64::Mul(Y[r][c - 1], p[r][YY]);
		}

	Fixed64::FP_LONG result = 0;
	for (int c = 1; c < COLS; c++) {
		result += M[ROWS - 1][c];
		result += X[ROWS - 1][c];
	}

	double result_double = log10(Fixed64::ToDouble(result)) - LOG10_INITIAL_CONSTANT;
	// std::cout << "Fixed64\t" << Fixed64::ToDouble(INITIAL_CONSTANT) << '\t' << LOG10_INITIAL_CONSTANT << '\t' << result_double << std::endl;;

	return result_double;
}
*/

void computeLikelihoodsNative_concurrent_trie(std::vector<trie_testcase> &testcases,
											  std::vector<std::vector<double>> &likelihoodArray) {
	for (int i = 0; i < testcases.size(); i++) {
		if (BOOST_LIKELY((!PairHMMConcurrentControl::startPairHMMConcurrentMode) || testcases.size() - i < 96)) {
			computeLikelihoodsNative_concurrent_trie_i(testcases, likelihoodArray, i);
		} else {
			std::shared_ptr<LikelihoodsTask_trie> likelihoods = std::make_shared<LikelihoodsTask_trie>(testcases, likelihoodArray, (unsigned long)i, testcases.size());
			PairHMMConcurrentControl::pairHMMMutex.lock();
			PairHMMConcurrentControl::pairHMMTaskQueue_trie.push(likelihoods);
			//std::cout << "push " + std::to_string(PairHMMConcurrentControl::pairHMMTaskQueue.size()) + '\n';
			PairHMMConcurrentControl::pairHMMMutex.unlock();
			//std::cout << std::to_string(testcases.size()) + " testcases\tpush [" + std::to_string(i) + ", " + std::to_string(testcases.size() - 1) + "] size: " +
			//                                                                                                                                          std::to_string(testcases.size() - i ) + '\n';

			for (unsigned long ind = likelihoods->index++; ind < likelihoods->testcasesSize - 8; ind = likelihoods->index++) {
				//std::cout << std::to_string(ind) + '\n';
				computeLikelihoodsNative_concurrent_trie_i(likelihoods->taskTestcases, likelihoods->taskLikelihoodArray, ind);
				likelihoods->count++;
			}

			for (unsigned long ind = likelihoods->testcasesSize - 8; ind < likelihoods->testcasesSize; ind++) {
				//std::cout << std::to_string(ind) + '\n';
				computeLikelihoodsNative_concurrent_trie_i(likelihoods->taskTestcases, likelihoods->taskLikelihoodArray, ind);
				likelihoods->count++;
			}

			// make sure all calculations have been completed
			while (likelihoods->count != likelihoods->testcasesSize)
				std::this_thread::yield();

			//std::cout << std::to_string(testcases.size()) + " testcases done.\n";
			return;
		}
	}
}

void computeLikelihoodsNative_concurrent_trie_i(std::vector<trie_testcase> &testcases,
                                                std::vector<std::vector<double>> &likelihoodArray, unsigned long i) {
	trie_testcase *cur_case = &testcases[i];
	std::vector<float> result_float = (BOOST_UNLIKELY(g_use_double)) ? std::vector<float>() : g_compute_full_prob_t_float(cur_case);
	for (int k = 0; k < result_float.size(); k++) {
		if (result_float[k] < MIN_ACCEPTED) {
            testcase ts(cur_case->haplotypeDataArray[k].length, cur_case->haplotypeDataArray[k].haplotypeBases, cur_case->readForPairHmm);
            double result = g_compute_full_prob_double(&ts);
            likelihoodArray[i].emplace_back(log10(result) - Context<double>::LOG10_INITIAL_CONSTANT);
		} else {
            likelihoodArray[i].emplace_back((double)(log10f(result_float[k]) - Context<float>::LOG10_INITIAL_CONSTANT));
        }
	}
}
