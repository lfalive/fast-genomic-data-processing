//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_SEQVERTEX_H
#define MUTECT2CPP_MASTER_SEQVERTEX_H

#include <utility>

#include "BaseVertex.h"

class SeqVertex : public BaseVertex {
public:
	/**
	 * Create a new SeqVertex with sequence and the next available id
	 * @param sequence our base sequence
	 */
	SeqVertex(const std::shared_ptr<uint8_t[]> &sequence, int length) : BaseVertex(sequence, length) {}

	SeqVertex() = default;

	[[nodiscard]] int getId() const { return (int) hashCode(); }

	[[nodiscard]] long hashCode() const;

	bool operator<(const SeqVertex &other) const;

	std::shared_ptr<SeqVertex> withoutSuffix(const std::shared_ptr<uint8_t[]> &suffix, int length);

	std::shared_ptr<SeqVertex>
	withoutPrefixAndSuffix(const std::shared_ptr<uint8_t[]> &prefix, int preLength,
	                       const std::shared_ptr<uint8_t[]> &suffix,
	                       int sufLength);
};


#endif //MUTECT2CPP_MASTER_SEQVERTEX_H
