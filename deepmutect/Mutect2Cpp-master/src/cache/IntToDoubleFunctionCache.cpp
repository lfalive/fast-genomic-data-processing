//
// Created by 梦想家xixi on 2021/10/30.
//

#include "IntToDoubleFunctionCache.h"
#include "Mutect2Utils.h"
#include <cmath>
#include <cstring>

double IntToDoubleFunctionCache::get(int i) {
	Mutect2Utils::validateArg(i >= 0, "Cache doesn't apply to negative number");
	if (i >= length) {
		if (i >= maxSize()) return compute(i);
		int newCapacity = std::max(i + 10, 2 * length);
		expandCache(newCapacity);
	}
	return cache[i];
}

//TODO:synchronized method
void IntToDoubleFunctionCache::expandCache(int newCapacity) {
	if (newCapacity < length)
		return;

	auto *newCache = new double[newCapacity + 1];
	if (length != 0)
		memcpy(newCache, cache, length * sizeof(double));

	for (int i = length; i < newCapacity + 1; i++) {
		newCache[i] = compute(i);
	}

	delete[] cache;
	cache = newCache;
	length = newCapacity + 1;
}

IntToDoubleFunctionCache::~IntToDoubleFunctionCache() {
	delete[] cache;
}
