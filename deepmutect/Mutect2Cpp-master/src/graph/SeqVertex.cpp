//
// Created by 梦想家xixi on 2021/11/15.
//

#include <cstring>
#include <memory>
#include "SeqVertex.h"
#include "graph/utils/GraphObjectPool.h"

long SeqVertex::hashCode() const {
	return (long) this;
}

bool SeqVertex::operator<(const SeqVertex &other) const {
	return hashCode() < other.hashCode();
}

std::shared_ptr<SeqVertex> SeqVertex::withoutSuffix(const std::shared_ptr<uint8_t[]> &suffix, int length) {
	int prefixSize = getLength() - length;
	int newLength = std::min(prefixSize, getLength());
	if (prefixSize > 0) {
		std::shared_ptr<uint8_t[]> newSequence(new uint8_t[newLength]);
		memcpy(newSequence.get(), sequence.get(), newLength);
		return GraphObjectPool::createSeqVertex(newSequence, newLength);
	}
	return nullptr;
}

std::shared_ptr<SeqVertex> SeqVertex::withoutPrefixAndSuffix(const std::shared_ptr<uint8_t[]> &prefix, int preLength,
                                                             const std::shared_ptr<uint8_t[]>& suffix, int sufLength) {
	int start = preLength;
	int length = getLength() - sufLength - preLength;
	int stop = start + length;
	int newLength = stop - start;
	if (length > 0) {
		std::shared_ptr<uint8_t[]> newSequence(new uint8_t[newLength]);
		memcpy(newSequence.get(), sequence.get() + start, newLength);
		return GraphObjectPool::createSeqVertex(newSequence, newLength);
	}
	return nullptr;
}
