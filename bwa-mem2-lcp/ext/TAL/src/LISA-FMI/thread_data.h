#ifndef THREAD_DATA_H
#define THREAD_DATA_H
#include"lisa_util.h"
class threadData {
	public:
		uint64_t *str_enc;
		int64_t *intv_all;
		int64_t numSMEMs;
		threadData(int64_t pool_size);
		void dealloc_td();
};

#endif
