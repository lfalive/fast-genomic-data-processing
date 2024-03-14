//
// Created by 梦想家xixi on 2021/11/9.
//

#include <cstring>
#include <utility>
#include <cassert>
#include "Haplotype.h"
#include "read/AlignmentUtils.h"
#include "SimpleInterval.h"

Haplotype::Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length, bool isRef) : Allele(copyArray(bases, length),
                                                                                               length, isRef),
                                                                                        eventMap(nullptr) {}

Haplotype::Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length, bool isRef, double score)
		: Allele(copyArray(bases, length), length, isRef), score(score), eventMap(nullptr) {}

std::shared_ptr<uint8_t[]> Haplotype::copyArray(const std::shared_ptr<uint8_t[]> &base, int length) {
	std::shared_ptr<uint8_t[]> res{new uint8_t[length + 1]{0}};
	memcpy(res.get(), base.get(), length);
	return res;
}

Haplotype::Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length) : Allele(copyArray(bases, length), length,
                                                                                   false), eventMap(nullptr) {}

void Haplotype::setCigar(std::shared_ptr<Cigar> &cigar) {
	this->cigar = AlignmentUtils::consolidateCigar(cigar);
	Mutect2Utils::validateArg(this->cigar->getReadLength() == getLength(),
	                          "Read length is not equal to the read length of the cigar");
}

Haplotype::Haplotype(const std::shared_ptr<uint8_t[]> &bases, bool isRef, int length, int alignmentStartHapwrtRef,
                     std::shared_ptr<Cigar> &cigar) : Allele(copyArray(bases, length), length, false),
                                                      alignmentStartHapwrtRef(alignmentStartHapwrtRef),
                                                      eventMap(nullptr) {
	setCigar(cigar);
}

Haplotype::Haplotype(const std::shared_ptr<uint8_t[]> &bases, int length, std::shared_ptr<Locatable> loc) : Allele(
		copyArray(bases, length), length, false), genomeLocation(std::move(loc)), eventMap(nullptr) {}

std::shared_ptr<Haplotype> Haplotype::trim(const std::shared_ptr<Locatable> &loc) {
	Mutect2Utils::validateArg(loc != nullptr, "Loc cannot be null");
	Mutect2Utils::validateArg(genomeLocation != nullptr, "Cannot trim a Haplotype without containing GenomeLoc");
	SimpleInterval interval = SimpleInterval(genomeLocation);
	Mutect2Utils::validateArg(interval.contains(loc), "Can only trim a Haplotype to a containing span.");
	Mutect2Utils::validateArg(getCigar() != nullptr, "Cannot trim haplotype without a cigar");

	int newStart = loc->getStart() - this->genomeLocation->getStart();
	int newStop = newStart + loc->getEnd() - loc->getStart();
	//std::cout << loc->getStart() + 1 << " " << loc->getEnd() + 1 << std::endl;
	std::pair<int, std::shared_ptr<uint8_t[]>> newBases = AlignmentUtils::getBasesCoveringRefInterval(newStart, newStop,
	                                                                                                  getBases(),
	                                                                                                  getBasesLength(),
	                                                                                                  0, getCigar());
	std::shared_ptr<Cigar> newCigar = AlignmentUtils::trimCigarByReference(getCigar(), newStart, newStop);

	if (newBases.second == nullptr || AlignmentUtils::startsOrEndsWithInsertionOrDeletion(newCigar))
		return nullptr;

	std::shared_ptr<Haplotype> ret = std::make_shared<Haplotype>(newBases.second, newBases.first, getIsReference());
	ret->setCigar(newCigar);
	ret->setGenomeLocation(loc);
	ret->setScore(score);
	ret->setAlignmentStartHapwrtRef(newStart + getAlignmentStartHapwrtRef());
	return ret;
}

std::shared_ptr<Cigar> Haplotype::getCigar() {
	return cigar;
}

void Haplotype::setGenomeLocation(const std::shared_ptr<Locatable> &genomeLocation) {
	this->genomeLocation = genomeLocation;
}

void Haplotype::setScore(double score) {
	this->score = score;
}

int Haplotype::getAlignmentStartHapwrtRef() const {
	return alignmentStartHapwrtRef;
}

EventMap *Haplotype::getEventMap() {
	return eventMap;
}

void Haplotype::setEventMap(EventMap *_eventMap) {
	this->eventMap = _eventMap;
}

bool Haplotype::operator<(const Haplotype &other) const {
	if (this->getLength() < other.getLength())
		return true;
	else if (this->getLength() > other.getLength())
		return false;
	else {
		uint8_t *bases = this->getBases().get();
		uint8_t *otherBases = other.getBases().get();
		for (int i = 0; i < this->getLength(); i++) {
			if (bases[i] < otherBases[i])
				return true;
			else if (bases[i] > otherBases[i])
				return false;
			else {
				continue;
			}
		}
		return false;
	}
}

double Haplotype::getScore() const {
	return score;
}

std::shared_ptr<Cigar> Haplotype::getConsolidatedPaddedCigar(int padSize) {
	assert(padSize >= 0);
	auto extendedHaplotypeCigar = std::make_shared<Cigar>(getCigar()->getCigarElements());
	if (padSize > 0) {
		extendedHaplotypeCigar->add(CigarElement(padSize, CigarOperator::M));
	}
	return AlignmentUtils::consolidateCigar(extendedHaplotypeCigar);
}
