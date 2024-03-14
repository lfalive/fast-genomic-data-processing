//
// Created by hlf on 8/16/22.
//

#ifndef MUTECT2CPP_MASTER_READFORPAIRHMM_H
#define MUTECT2CPP_MASTER_READFORPAIRHMM_H

#include "xxhash.hpp"
#include "parallel_hashmap/phmap.h"

struct ReadForPairHMM {
	int rslen;
	const uint8_t *rs;
	uint32_t *charCombination;
	xxh::hash64_t hashCode;

	std::shared_ptr<float> floatInitializedVectors;

	ReadForPairHMM(int _rslen, const uint8_t *readQuals, const uint8_t *insGops, const uint8_t *delGops,
	               const char *gapConts, const uint8_t *reads);

	~ReadForPairHMM();

	template<typename NUMBER>
	std::shared_ptr<NUMBER> initializeData();

	template<typename NUMBER>
	std::shared_ptr<NUMBER> getInitializedData();

	void initializeFloatVector();
};

struct ReadForPairHMMHash {
	xxh::hash64_t operator()(const std::shared_ptr<ReadForPairHMM> &t) const {
		return t->hashCode;
	}
};

struct ReadForPairHMMEqual {
	bool operator()(const std::shared_ptr<ReadForPairHMM> &t1, const std::shared_ptr<ReadForPairHMM> &t2) const {
		if (t1->rslen != t2->rslen)
			return false;

		if (memcmp(t1->rs, t2->rs, t1->rslen) != 0)
			return false;

		return memcmp(t1->charCombination, t2->charCombination, 4 * t1->rslen) == 0;
	}
};


#endif //MUTECT2CPP_MASTER_READFORPAIRHMM_H
