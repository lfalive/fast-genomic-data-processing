//
// Created by 梦想家xixi on 2021/10/19.
//

#ifndef MUTECT2CPP_MASTER_KMER_H
#define MUTECT2CPP_MASTER_KMER_H


#include <cstdint>
#include <memory>
#include "xxhash.hpp"

class Kmer {
private:

	std::shared_ptr<uint8_t[]> bases;
	int start;
	int length;
	xxh::hash64_t hash;

	/**
	 *  Compute the hashcode for a KMer.
	 *  Equivalent to <code>new String(bases, start, length).hashCode()</code>
	 */
	static xxh::hash64_t hashCode(const std::shared_ptr<uint8_t[]> &bases, int start, int length);

public:
	Kmer(std::shared_ptr<uint8_t[]> kmer, int length);

	Kmer(Kmer const &kmer);

	/**
	 * Create a new kmer backed by the bases in bases, spanning start -> start + length
	 *
	 * Under no circumstances can bases be modified anywhere in the client code.  This does not make a copy
	 * of bases for performance reasons
	 *
	 * @param bases an array of bases
	 * @param start the start of the kmer in bases, must be >= 0 and < bases.length
	 * @param length the length of the kmer.  Must be >= 0 and start + length < bases.length
	 */
	Kmer(std::shared_ptr<uint8_t[]> kmer, int start, int length);

	Kmer subKmer(int newStart, int newLength);

	//delete
	[[nodiscard]] std::shared_ptr<uint8_t[]> getBases() const;

	/**
	 * The length of this kmer
	 * @return an integer >= 0
	 */
	[[nodiscard]] int getLength() const { return length; }

	int getDifferingPositions(const Kmer &other, int maxDistance, const std::shared_ptr<int> &differingIndeces,
	                          const std::shared_ptr<uint8_t[]> &differingBases);

	bool operator<(const Kmer &other) const;

	bool operator==(const Kmer &other) const;

	[[nodiscard]] xxh::hash64_t getHash() const { return hash; };

};

struct equal_kmer {
	bool operator()(const std::shared_ptr<Kmer> &kmer1, const std::shared_ptr<Kmer> &kmer2) const;
};

struct hash_kmer {
	xxh::hash64_t operator()(const std::shared_ptr<Kmer> &kmer1) const;
};


#endif //MUTECT2CPP_MASTER_KMER_H
