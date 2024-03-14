//
// Created by 梦想家xixi on 2021/10/19.
//

#include "Kmer.h"
#include "Mutect2Utils.h"
#include <cstring>
#include <utility>

Kmer::Kmer(std::shared_ptr<uint8_t[]> kmer, const int length) : bases(std::move(kmer)), start(0), length(length) {
	Mutect2Utils::validateArg(start >= 0, "start must be >= 0");
	Mutect2Utils::validateArg(length >= 0, "length must be >= 0");
	this->hash = hashCode(bases, start, length);
}

Kmer::Kmer(std::shared_ptr<uint8_t[]> kmer, const int start, const int length) : bases(std::move(kmer)), start(start),
                                                                                 length(length) {
	Mutect2Utils::validateArg(start >= 0, "start must be >= 0");
	Mutect2Utils::validateArg(length >= 0, "length must be >= 0");
	this->hash = hashCode(bases, start, length);
}

Kmer::Kmer(const Kmer &kmer) : bases(kmer.bases), start(kmer.start), length(kmer.length) {
	Mutect2Utils::validateArg(start >= 0, "start must be >= 0");
	Mutect2Utils::validateArg(length >= 0, "length must be >= 0");
	this->hash = hashCode(bases, start, length);
}

xxh::hash64_t Kmer::hashCode(const std::shared_ptr<uint8_t[]> &bases, const int start, const int length) {
	if (length == 0) {
		return 0;
	}
	//int h = 0;
	uint8_t *bases_ = bases.get();
	return xxh::xxhash3<64>(bases_ + start, length);
	/* for(int i = start, stop = start + length; i < stop; i++) {
		 h = 31 * h + bases_[i];
	 }
	 return h;*/
}

Kmer Kmer::subKmer(const int newStart, const int newLength) {
	Kmer subkmer = Kmer(bases, start + newStart, newLength);
	return subkmer;
}

std::shared_ptr<uint8_t[]> Kmer::getBases() const {
	if (length == 0)
		return nullptr;
	std::shared_ptr<uint8_t[]> res(new uint8_t[length]);
	memcpy(res.get(), bases.get() + start, sizeof(uint8_t) * length);
	return res;
}

int Kmer::getDifferingPositions(const Kmer &other, int maxDistance, const std::shared_ptr<int> &differingIndeces,
                                const std::shared_ptr<uint8_t[]> &differingBases) {
	Mutect2Utils::validateArg(differingIndeces != nullptr, "Null object is not allowed here.");
	Mutect2Utils::validateArg(differingBases != nullptr, "Null object is not allowed here.");
	Mutect2Utils::validateArg(maxDistance > 0, "maxDistance must be positive");
	int dist = 0;
	if (length == other.length) {
		uint8_t *f2 = other.getBases().get();
		uint8_t *bases_ = bases.get();
		int *differingIndeces_ = differingIndeces.get();
		uint8_t *differingBases_ = differingBases.get();
		for (int i = 0; i < length; i++) {
			if (bases_[start + i] != f2[i]) {
				differingIndeces_[dist] = i;
				differingBases_[dist++] = f2[i];
				if (dist > maxDistance) {
					return -1;
				}
			}
		}
	}
	return dist;
}

bool Kmer::operator<(const Kmer &other) const {
	if (this->length > other.length)
		return false;
	if (this->length < other.length)
		return true;
	return memcmp(bases.get() + this->start, other.bases.get() + other.start, length) < 0;
}

bool Kmer::operator==(const Kmer &other) const {
	if (this->hash != other.hash || this->length != other.length)
		return false;
	return memcmp(bases.get() + this->start, other.bases.get() + other.start, length) == 0;
}

bool equal_kmer::operator()(const std::shared_ptr<Kmer> &kmer1, const std::shared_ptr<Kmer> &kmer2) const {
	return *kmer1 == *kmer2;
}

xxh::hash64_t hash_kmer::operator()(const std::shared_ptr<Kmer> &kmer1) const {
	return kmer1->getHash();
}
