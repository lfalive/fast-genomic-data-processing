//
// Created by 梦想家xixi on 2021/11/9.
//

#ifndef MUTECT2CPP_MASTER_HAPLOTYPE_H
#define MUTECT2CPP_MASTER_HAPLOTYPE_H


#include "Allele.h"
#include "Locatable.h"
#include "cigar/Cigar.h"
#include "EventMap.h"


class EventMap;

class Haplotype : public Allele {
private:
	std::shared_ptr<Locatable> genomeLocation;
	std::shared_ptr<Cigar> cigar;
	EventMap *eventMap;
	int alignmentStartHapwrtRef{};
	double score{};

	static std::shared_ptr<uint8_t[]> copyArray(const std::shared_ptr<uint8_t[]> &base, int length);

public:
	/**
	 * Main constructor
	 *
	 * @param bases a non-null array of bases
	 * @param isRef is this the reference haplotype?
	 */
	Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length, bool isRef);

	/**
	 * Create a new non-ref haplotype
	 *
	 * @param bases a non-null array of bases
	 */
	Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length);

	/**
	* Create a new haplotype with bases
	*
	* Requires bases.length == cigar.getReadLength()
	*
	* @param bases a non-null array of bases
	* @param isRef is this the reference haplotype?
	* @param alignmentStartHapwrtRef offset of this haplotype w.r.t. the reference
	* @param cigar the cigar that maps this haplotype to the reference sequence
	*/
	Haplotype(const std::shared_ptr<uint8_t[]> &bases, bool isRef, int length, int alignmentStartHapwrtRef,
	          std::shared_ptr<Cigar> &cigar);

	Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length, std::shared_ptr<Locatable> loc);

	Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length, bool isRef, double score);

	~Haplotype() = default;

	/**
	* Set the cigar of this haplotype to cigar.
	*
	* Note that this function consolidates the cigar, so that 1M1M1I1M1M => 2M1I2M
	*
	* @param cigar a cigar whose readLength == length()
	*/
	void setCigar(std::shared_ptr<Cigar> &cigar);

	/**
	* Create a new Haplotype derived from this one that exactly spans the provided location
	*
	* Note that this haplotype must have a contain a genome loc for this operation to be successful.  If no
	* GenomeLoc is contained than @throws an IllegalStateException
	*
	* Also loc must be fully contained within this Haplotype's genomeLoc.  If not an IllegalArgumentException is
	* thrown.
	*
	* @param loc a location completely contained within this Haplotype's location
	* @return a new Haplotype within only the bases spanning the provided location, or null for some reason the haplotype would be malformed if
	*/
	std::shared_ptr<Haplotype> trim(const std::shared_ptr<Locatable> &loc);

	/**
	* Get the cigar for this haplotype.  Note that the cigar is guaranteed to be consolidated
	* in that multiple adjacent equal operates will have been merged
	* @return the cigar of this haplotype
	*/
	std::shared_ptr<Cigar> getCigar();

	void setGenomeLocation(const std::shared_ptr<Locatable> &genomeLocation);

	void setScore(double score);

	void setAlignmentStartHapwrtRef(int alignmentStartHapwrtRef) {
		this->alignmentStartHapwrtRef = alignmentStartHapwrtRef;
	}

	int getAlignmentStartHapwrtRef() const;

	const std::shared_ptr<Locatable> &getGenomeLocation() { return genomeLocation; }

	EventMap *getEventMap();

	void setEventMap(EventMap *_eventMap);

	bool operator<(const Haplotype &other) const;

	double getScore() const;


	/**
    * Get the haplotype cigar extended by padSize M at the tail, consolidated into a clean cigar
    *
    * @param padSize how many additional Ms should be appended to the end of this cigar.  Must be >= 0
    * @return a newly allocated Cigar that consolidate(getCigar + padSize + M)
    */
	std::shared_ptr<Cigar> getConsolidatedPaddedCigar(int padSize);
};


#endif //MUTECT2CPP_MASTER_HAPLOTYPE_H
