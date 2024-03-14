//
// Created by 梦想家xixi on 2021/10/11.
//

#ifndef MUTECT2CPP_MASTER_SIMPLEINTERVAL_H
#define MUTECT2CPP_MASTER_SIMPLEINTERVAL_H

#include "Mutect2Utils.h"
#include "IntervalUtils.h"
#include <string>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include "samtools/SAMSequenceDictionary.h"
#include "Locatable.h"
#include "ContigMap.h"

static const char CONTIG_SEPARATOR = ':';
static const char START_END_SEPARATOR = '-';
static const std::string END_OF_CONTIG = "+";

class SimpleInterval : public Locatable {
private:
	int start;
	int end;
	int contig;

	bool contiguous(Locatable *other) const;

public:
	SimpleInterval(const std::string &contig, int start, int end);

	SimpleInterval(std::string &&contig, int start, int end);

	SimpleInterval(int contig, int start, int end);

	SimpleInterval(SimpleInterval const &simpleInterval);

	/**
	 * Makes an interval by parsing the string.
	 *
	 * @warning this method does not fill in the true contig end values
	 * for intervals that reach to the end of their contig,
	 * uses {@link Integer#MAX_VALUE} instead.
	 *
	 * Semantics of start and end are defined in {@link Locatable}.
	 * The format is one of:
	 *
	 * contig           (Represents the whole contig, from position 1 to the {@link Integer#MAX_VALUE})
	 *
	 * contig:start     (Represents the 1-element range start-start on the given contig)
	 *
	 * contig:start-end (Represents the range start-end on the given contig)
	 *
	 * contig:start+    (Represents the prefix of the contig starting at the given start position and ending at {@link Integer#MAX_VALUE})
	 *
	 * examples (note that _all_ commas in numbers are simply ignored, for human convenience):
	 *
	 * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000' or 'chr2:1000000+'
	  *
	  * @param str non-empty string to be parsed
	 */
	SimpleInterval(std::string &str);

	SimpleInterval();

	virtual ~SimpleInterval() = default;

	void clearContig();

	/**
	 * Parses a number like 100000 or 1,000,000 into an int.
	 */
	static int parsePosition(std::string pos);

	bool operator==(const SimpleInterval &interval) const;

	bool equal(const SimpleInterval &interval) const { return *this == interval; }

	std::string getContig() const override { return ContigMap::getContigString(contig); }

	int getContigInt() const { return contig; }

	int getStart() const override { return start; }

	int getEnd() const override { return end; }

	int size() const { return end - start + 1; }

	/**
	  * Determines whether this interval comes within "margin" of overlapping the provided locatable.
	  * This is the same as plain overlaps if margin=0.
	  *
	  * @param other interval to check
	  * @param margin how many bases may be between the two interval for us to still consider them overlapping; must be non-negative
	  * @return true if this interval overlaps other, otherwise false
	  * @throws IllegalArgumentException if margin is negative
	  */
	bool overlapsWithMargin(const std::shared_ptr<SimpleInterval> &other, int margin) const;

	/**
	 * Determines whether this interval overlaps the provided locatable.
	 *
	 * @param other interval to check
	 * @return true if this interval overlaps other, otherwise false
	 */
	bool overlaps(const std::shared_ptr<SimpleInterval> &other) const;

	/**
	  * Returns the intersection of the two intervals. The intervals must overlap or IllegalArgumentException will be thrown.
	 */
	std::shared_ptr<SimpleInterval> intersect(const std::shared_ptr<SimpleInterval> &other) const;

	/**
	  * Returns a new SimpleInterval that represents the entire span of this and other.  Requires that
	  * this and that SimpleInterval are contiguous.
	  */
	std::shared_ptr<SimpleInterval> mergeWithContiguous(const std::shared_ptr<SimpleInterval> &other);

	/**
	  * Returns a new SimpleInterval that represents the region between the endpoints of this and other.
	  *
	  * Unlike {@link #mergeWithContiguous}, the two intervals do not need to be contiguous
	  *
	  * @param other the other interval with which to calculate the span
	  * @return a new SimpleInterval that represents the region between the endpoints of this and other.
	  */
	std::shared_ptr<SimpleInterval> spanWith(const std::shared_ptr<SimpleInterval> &other) const;

	/**
	  * Returns a new SimpleInterval that represents this interval as expanded by the specified amount in both
	  * directions, bounded by the contig start/stop if necessary.
	  *
	  * @param padding amount to expand this interval
	  * @param contigLength length of this interval's contig
	  * @return a new SimpleInterval that represents this interval as expanded by the specified amount in both
	  *         directions, bounded by the contig start/stop if necessary.
	  */
	std::shared_ptr<SimpleInterval> expandWithinContig(int padding, int contigLength) const;

	std::shared_ptr<SimpleInterval> expandWithinContig(int padding, SAMSequenceDictionary *sequenceDictionary) const;

	friend std::ostream &operator<<(std::ostream &os, const SimpleInterval &simpleInterval);

	SimpleInterval(const std::shared_ptr<Locatable> &pLocatable);

	void printInfo() const;

	bool contiguous(SimpleInterval *other) const;
};

#endif //MUTECT2CPP_MASTER_SIMPLEINTERVAL_H