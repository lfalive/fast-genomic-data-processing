//
// Created by 梦想家xixi on 2021/10/20.
//

#include "BaseVertex.h"
#include "Mutect2Utils.h"

xxh::hash64_t BaseVertex::hashCode(const std::shared_ptr<uint8_t[]> &a, int length) {
	if (a == nullptr)
		return 0;
	return xxh::xxhash3<64>(a.get(), length);
}

BaseVertex::BaseVertex(std::shared_ptr<uint8_t[]> const &sequence, const int length) : sequence(sequence),
                                                                                       length(length) {
	Mutect2Utils::validateArg(sequence != nullptr || length == 0, "Sequence cannot be null");
	cashedHashCode = hashCode(sequence, length);
}

bool BaseVertex::isEmpty() const {
	return length == 0;
}

bool BaseVertex::operator==(const BaseVertex &other) const {
	if (other.cashedHashCode != cashedHashCode || other.length != length)
		return false;
	return memcmp(sequence.get(), other.sequence.get(), length) == 0;
}

bool BaseVertex::operator<(const BaseVertex &other) const {
	if (length >= other.length)
		return false;
	for (int i = 0; i < length; i++)
		if (sequence.get()[i] > other.sequence.get()[i])
			return false;
	return true;
}

std::ostream &operator<<(std::ostream &os, const BaseVertex &baseVertex) {
	os << "baseVertex : ";
	for (int i = 0; i < baseVertex.length; i++)
		os << baseVertex.sequence.get()[i];
	os << '.' << std::endl;
	return os;
}

void BaseVertex::setAdditionalInfo(const std::string &info) {
	additionalInfo = info;
}

void BaseVertex::additionalInfoAppendPlusSign() {
	additionalInfo += plusSign;
}


bool BaseVertex::hasAmbiguousSequence() {
	for (int i = 0; i < length; i++) {
		uint8_t tmp = sequence.get()[i];
		if (tmp > 60)
			tmp -= 32;
		switch (tmp) {
			case 'A':
			case 'T':
			case 'G':
			case 'C':
				continue;
			default:
				return true;
		}
	}
	return false;
}

bool BaseVertex::seqEquals(const std::shared_ptr<BaseVertex> &other) {
	if (length != other->getLength() || cashedHashCode != other->getHashCode())
		return false;

	return memcmp(sequence.get(), other->getSequence().get(), length) == 0;
}
