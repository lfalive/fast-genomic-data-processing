#include "thread_data.h"

threadData::threadData(int64_t pool_size){

	numSMEMs = 0;

	str_enc = (uint64_t *)aligned_alloc(64, pool_size * sizeof(uint64_t));
	intv_all = (int64_t *)aligned_alloc(64, pool_size * 2 * sizeof(int64_t));
}

void threadData::dealloc_td(){
	free(str_enc);
	free(intv_all);
}

