//
// Created by 梦想家xixi on 2021/10/20.
//

#ifndef MUTECT2CPP_MASTER_BASEVERTEX_H
#define MUTECT2CPP_MASTER_BASEVERTEX_H


#include <cstdint>
#include <string>
#include <iostream>
#include <memory>
#include "xxhash.hpp"

static char plusSign = '+';

class BaseVertex {
private:
	xxh::hash64_t cashedHashCode;
	int length;
	std::string additionalInfo;

	static xxh::hash64_t hashCode(const std::shared_ptr<uint8_t[]> &a, int length);

protected:
	std::shared_ptr<uint8_t[]> sequence;

public:
	/**
	 * Create a new sequence vertex with sequence
	 *
	 * This code doesn't copy sequence for efficiency reasons, so sequence must absolutely not be modified
	 * in any way after passing this sequence to the BaseVertex
	 *
	 * @param sequence a non-null sequence of bases contained in this vertex
	 */
	BaseVertex(const std::shared_ptr<uint8_t[]> &sequence, int length);

	BaseVertex() = default;

	virtual ~BaseVertex() = default;

	/**
	 * Does this vertex have an empty sequence?
	 *
	 * That is, is it a dummy node that's only present for structural reasons but doesn't actually
	 * contribute to the sequence of the graph?
	 *
	 * @return true if sequence is empty, false otherwise
	 */
	[[nodiscard]] bool isEmpty() const;

	/**
	 * Get the length of this sequence
	 * @return a positive integer >= 1
	 */
	[[nodiscard]] int getLength() const { return length; }

	bool operator==(const BaseVertex &other) const;

	bool operator<(const BaseVertex &other) const;

	[[nodiscard]] xxh::hash64_t getHashCode() const { return cashedHashCode; }

	friend std::ostream &operator<<(std::ostream &os, const BaseVertex &baseVertex);

	[[nodiscard]] std::shared_ptr<uint8_t[]> getSequence() const { return sequence; }

	void setAdditionalInfo(const std::string &info);

	[[nodiscard]] std::string getAdditionalInfo() const { return additionalInfo; }

	bool hasAmbiguousSequence();

	bool seqEquals(const std::shared_ptr<BaseVertex> &other);

	virtual std::shared_ptr<uint8_t[]> getAdditionalSequence(bool source) { return getSequence(); }

	virtual int getAdditionalSequenceLength(bool source) { return getLength(); }

	void additionalInfoAppendPlusSign();
};


#endif //MUTECT2CPP_MASTER_BASEVERTEX_H
