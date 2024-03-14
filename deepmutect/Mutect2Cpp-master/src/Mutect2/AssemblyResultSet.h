//
// Created by 梦想家xixi on 2021/12/1.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
#define MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H

#include "parallel_hashmap/phmap.h"
#include "AssemblyResult.h"
#include "Haplotype.h"
#include "AssemblyRegion.h"
#include "VariantContext.h"

struct HaplotypeComp {
public:
	bool operator()(const std::shared_ptr<Haplotype> &left, const std::shared_ptr<Haplotype> &right) const {
		if (left->getLength() != right->getLength())
			return left->getLength() > right->getLength();
		return memcmp(left->getBases().get(), right->getBases().get(), left->getLength()) < 0;
	}
};

struct hash_Haplotype {
	size_t operator()(const std::shared_ptr<Haplotype> &haplotype) const {
		return xxh::xxhash3<64>(haplotype->getBases().get(), haplotype->getBasesLength());
	}
};

struct equal_Haplotype {
	bool operator()(const std::shared_ptr<Haplotype> &left, const std::shared_ptr<Haplotype> &right) const {
		if (left->getLength() != right->getLength())
			return false;
		return memcmp(left->getBases().get(), right->getBases().get(), left->getLength()) == 0;
	}
};

class AssemblyResultSet {
private:
	std::map<int, std::shared_ptr<AssemblyResult>> assemblyResultByKmerSize;
	std::set<std::shared_ptr<Haplotype>, HaplotypeComp> haplotypes;
	std::map<std::shared_ptr<Haplotype>, std::shared_ptr<AssemblyResult>, HaplotypeComp> assemblyResultByHaplotype;
	std::shared_ptr<AssemblyRegion> regionForGenotyping;
	std::shared_ptr<uint8_t[]> fullReferenceWithPadding;
	int fullReferenceWithPaddingLength{};
	std::shared_ptr<SimpleInterval> paddedReferenceLoc;
	bool variationPresent{};
	std::shared_ptr<Haplotype> refHaplotype;
	bool wasTrimmed = false;
	int lastMaxMnpDistanceUsed = -1;
	std::set<int> kmerSizes;
	std::set<std::shared_ptr<VariantContext>, VariantContextComparator> variationEvents;

	bool add(const std::shared_ptr<AssemblyResult> &ar);

	void updateReferenceHaplotype(const std::shared_ptr<Haplotype> &newHaplotype);

	std::vector<std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>>> calculateOriginalByTrimmedHaplotypes(
			const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion);

	static std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
	trimDownHaplotypes(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion,
	                   const std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> &haplotypeList);

public:
	AssemblyResultSet() {};

	bool add(const std::shared_ptr<Haplotype> &h, const std::shared_ptr<AssemblyResult> &ar);

	bool add(const std::shared_ptr<Haplotype> &h);

	void setRegionForGenotyping(std::shared_ptr<AssemblyRegion> regionForGenotyping);

	void setFullReferenceWithPadding(std::shared_ptr<uint8_t[]> fullReferenceWithPadding, int length);

	void setPaddedReferenceLoc(const std::shared_ptr<SimpleInterval> &paddedReferenceLoc);

	std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &getVariationEvents(int maxMnpDistance);

	void regenerateVariationEvents(int distance);

	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> getHaplotypeList();

	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> getSortedHaplotypeList();

	bool isisVariationPresent();

	void deleteEventMap();

	~AssemblyResultSet() = default;

	/**
	 * Trims an assembly result set down based on a new set of trimmed haplotypes.
	 *
	 * @param trimmedAssemblyRegion the trimmed down active region.
	 *
	 * @throws NullPointerException if any argument in {@code null} or
	 *      if there are {@code null} entries in {@code originalByTrimmedHaplotypes} for trimmed haplotype keys.
	 * @throws IllegalArgumentException if there is no reference haplotype amongst the trimmed ones.
	 *
	 * @return never {@code null}, a new trimmed assembly result set.
	 */
	std::shared_ptr<AssemblyResultSet> trimTo(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion);

	/**
	 * Returns the current region for genotyping.
	 *
	 * @return might be {@code null}.
	 */
	std::shared_ptr<AssemblyRegion> getRegionForGenotyping();

	void printSortedHaplotypes();

	/**
    * Query the reference haplotype in the result set.
    * @return {@code null} if none wasn't yet added, otherwise a reference haplotype.
    */
	std::shared_ptr<Haplotype>& getReferenceHaplotype();

	/**
    * Returns the padded reference location.
    *
    * @return might be {@code null}
    */
	std::shared_ptr<SimpleInterval>& getPaddedReferenceLoc();

	/**
    * Returns the current full reference with padding.
    *
    * @return might be {@code null}. The result must not be modified by the caller.
    */
	std::shared_ptr<uint8_t[]>& getFullReferenceWithPadding();

	int getFullReferenceWithPaddingLength();
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYRESULTSET_H
