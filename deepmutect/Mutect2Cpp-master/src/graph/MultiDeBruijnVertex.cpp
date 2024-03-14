//
// Created by 梦想家xixi on 2021/10/20.
//

#include "MultiDeBruijnVertex.h"

MultiDeBruijnVertex::MultiDeBruijnVertex(const std::shared_ptr<uint8_t[]> &sequence, int length,
                                         bool mergeIdenticalNodes)
		: BaseVertex(sequence, length), mergeIdenticalNodes(mergeIdenticalNodes) {
	hashCode = mergeIdenticalNodes ? (long) BaseVertex::getHashCode() : (long) this;
}

bool MultiDeBruijnVertex::operator==(const MultiDeBruijnVertex &other) const {
	if (this->getLength() != other.getLength() || this->getHashCode() != other.getHashCode())
		return false;
	return memcmp(sequence.get(), other.sequence.get(), this->getLength()) == 0;
}

bool MultiDeBruijnVertex::operator<(const MultiDeBruijnVertex &other) const {
	if (this->getLength() > other.getLength())
		return false;
	if (this->getLength() < other.getLength())
		return true;
	return memcmp(sequence.get(), other.sequence.get(), this->getLength()) <= 0;
}

std::shared_ptr<uint8_t[]> MultiDeBruijnVertex::getAdditionalSequence(bool source) {
	return source ? BaseVertex::getAdditionalSequence(source) : getSuffixAsArray();
}

std::shared_ptr<uint8_t[]> MultiDeBruijnVertex::getSuffixAsArray() const {
	std::shared_ptr<uint8_t[]> res(new uint8_t[1]);
	res.get()[0] = getSuffix();
	return res;
}

MultiDeBruijnVertex::MultiDeBruijnVertex(const std::shared_ptr<uint8_t[]> &sequence, int length)
		: BaseVertex(sequence, length), mergeIdenticalNodes(false) {
	hashCode = mergeIdenticalNodes ? (long) BaseVertex::getHashCode() : (long) this;
}

int MultiDeBruijnVertex::getAdditionalLength(bool source) {
	return source ? getLength() : 1;
}

int MultiDeBruijnVertex::getAdditionalSequenceLength(bool isSource) {
	return isSource ? getLength() : 1;
}
