//
// Created by hlf on 8/27/22.
//

#include <variant>
#include "ReadForPairHMM.h"
#include "intel/pairhmm/Context.h"
#include "intel/common/avx.h"


/**
 * In order to assign the memory address of the pre-processing array to the SIMD type pointer,
 * it is necessary to align the memory according to the SIMD_TYPE data.
 * */
static bool use_avx512 = is_avx512_supported();
static size_t SIMD_size = use_avx512 ? sizeof(__m512) : sizeof(__m256);

ReadForPairHMM::ReadForPairHMM(int _rslen, const uint8_t *readQuals, const uint8_t *insGops, const uint8_t *delGops,
                               const char *gapConts, const uint8_t *reads) : rslen(_rslen), rs(reads) {
	/**
	 * Calculate the hashcode and store it
	 * The calculation of hashcode does not consider the base sequence,
	 * which can remove most of the duplication and improve the efficiency.
	 * But the base sequence should be considered when determining whether the two are equal.
	 */

	// 4 * sizeof(uint8_t) * rslen = sizeof(uint32_t) * rslen
	charCombination = new uint32_t[rslen];
	auto *charP = (char *) charCombination;
	memcpy(charP, delGops, rslen);
	memcpy(charP + rslen, insGops, rslen);
	memcpy(charP + 2 * rslen, gapConts, rslen);
	memcpy(charP + 3 * rslen, readQuals, rslen);
	for (int k = 0; k < rslen; ++k) {
		charCombination[k] &= 2139062143;   // means every char &= 127
	}
	hashCode = xxh::xxhash3<64>(charCombination, rslen);
}

template<typename NUMBER>
std::shared_ptr<NUMBER> ReadForPairHMM::initializeData() {
	Context<NUMBER> ctx;
	int AVX_LENGTH = SIMD_size / sizeof(NUMBER);
	int SIMD_COUNT = (rslen + AVX_LENGTH) / AVX_LENGTH;
	int NUMBER_COUNT = SIMD_COUNT * AVX_LENGTH;

	/**
	 * In order to better adapt to SIMD registers, the actual array size is larger than ROW.
	 * For calculation of a certain precision, AVX_LENGTH is fixed.
	 * */

	std::shared_ptr<NUMBER> ret(
			use_avx512 ? (NUMBER *) new __m512[7 * SIMD_COUNT] : (NUMBER *) new __m256[7 * SIMD_COUNT],
			[](NUMBER *p) {
				if (use_avx512)
					delete[] (__m512*) p;
				else
					delete[] (__m256*) p;
			});

	NUMBER *p_MM = ret.get();
	NUMBER *p_XX = ret.get() + NUMBER_COUNT;
	NUMBER *p_YY = ret.get() + 2 * NUMBER_COUNT;
	NUMBER *p_MX = ret.get() + 3 * NUMBER_COUNT;
	NUMBER *p_MY = ret.get() + 4 * NUMBER_COUNT;
	NUMBER *p_GAPM = ret.get() + 5 * NUMBER_COUNT;
	NUMBER *distm1D = ret.get() + 6 * NUMBER_COUNT;
	auto *_d = (uint8_t *) charCombination;
	auto *_i = (uint8_t *) charCombination + rslen;
	auto *_c = (uint8_t *) charCombination + 2 * rslen;
	auto *_q = (uint8_t *) charCombination + 3 * rslen;
	for (int r = 0; r < rslen; r++) {
		*p_MM++ = ctx.set_mm_prob(*_i, *_d); // p_MM
		*p_XX++ = ctx.ph2pr[*_c];   // p_XX
		*p_YY++ = ctx.ph2pr[*_c];   // p_YY
		*p_MX++ = ctx.ph2pr[*_i++];   // p_MX
		*p_MY++ = ctx.ph2pr[*_d++];   // p_MY
		*p_GAPM++ = ctx._(1.0) - ctx.ph2pr[*_c++];  // p_GAPM
		*distm1D++ = ctx.ph2pr[*_q++];   // distm1D
	}
	return ret;
}

/**
 * In order to balance the cost of calculation and memory application:
 * For float calculation, store the result of pre calculation.
 * For double calculation, perform pre calculation and return the result during compute_full_prob.
 * */
template<>
std::shared_ptr<float> ReadForPairHMM::getInitializedData() {
	return floatInitializedVectors;
}

template<>
std::shared_ptr<double> ReadForPairHMM::getInitializedData() {
	return initializeData<double>();
}

void ReadForPairHMM::initializeFloatVector() {
	floatInitializedVectors = initializeData<float>();
}

ReadForPairHMM::~ReadForPairHMM() {
	delete[] charCombination;
}
