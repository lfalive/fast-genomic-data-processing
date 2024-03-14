//
// Created by 梦想家xixi on 2021/12/4.
//

#include "EventMap.h"

#include <utility>
#include <vector>
#include "utils/BaseUtils.h"
#include "variantcontext/builder/VariantContextBuilder.h"
#include <deque>
#include <iostream>

const std::shared_ptr<Allele> EventMap::SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele::create(std::shared_ptr<uint8_t[]>(
		                                                                                           new uint8_t[19]{'<', 'U', 'N', 'A', 'S', 'S', 'E', 'M', 'B', 'L', 'E', 'D', '_', 'E', 'V', 'E', 'N', 'T', '>'}),
                                                                                           19, false);

EventMap::EventMap(std::shared_ptr<Haplotype> haplotype, std::shared_ptr<uint8_t[]> ref, int refLength,
                   std::shared_ptr<Locatable> refLoc, std::string sourceNameToAdd,
                   int maxMnpDistance) : haplotype(std::move(haplotype)), ref(std::move(ref)), refLength(refLength),
                                         refLoc(std::move(refLoc)), sourceNameToAdd(std::move(sourceNameToAdd)) {
	processCigarForInitialEvents(maxMnpDistance);
}

void EventMap::processCigarForInitialEvents(int maxMnpDistance) {
	if (maxMnpDistance < 0)
		throw std::invalid_argument("maxMnpDistance may not be negative.");
	std::shared_ptr<Cigar> cigar = haplotype->getCigar();
	std::shared_ptr<uint8_t[]> alignment = haplotype->getBases();
	int alignmentLength = haplotype->getLength();

	int refPos = haplotype->getAlignmentStartHapwrtRef();
	if (refPos < 0) {
		return;
	}
	std::vector<std::shared_ptr<VariantContext>> proposedEvents;
	int alignmentPos = 0;
	for (int cigarIndex = 0; cigarIndex < cigar->numCigarElements(); cigarIndex++) {
		CigarElement ce = cigar->getCigarElement(cigarIndex);
		int elementLength = ce.getLength();
		switch (ce.getOperator()) {
			case I: {
				if (refPos > 0) {
					std::shared_ptr<std::vector<std::shared_ptr<Allele>>> insertionAlleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
					int insertionStart = refLoc->getStart() + refPos - 1;
					uint8_t refByte = ref.get()[refPos - 1];
					if (BaseUtils::isRegularBase(refByte)) {
						insertionAlleles->emplace_back(Allele::create(refByte, true));
					}
					if (cigarIndex == 0 || cigarIndex == cigar->numCigarElements() - 1) {

					} else {
						std::shared_ptr<uint8_t[]> insertionBases{new uint8_t[1]};
						insertionBases.get()[0] = ref.get()[refPos - 1];
						int toAddLength;
						std::shared_ptr<uint8_t[]> toAdd = Mutect2Utils::copyOfRange(alignment, alignmentLength,
						                                                             alignmentPos,
						                                                             alignmentPos + elementLength,
						                                                             toAddLength);
						int length = 1 + toAddLength;
						std::shared_ptr<uint8_t[]> tmp{new uint8_t[length]};
						memcpy(tmp.get(), insertionBases.get(), 1);
						memcpy(tmp.get() + 1, toAdd.get(), toAddLength);
						insertionBases = tmp;
						if (BaseUtils::isAllRegularBases(insertionBases, length)) {
							insertionAlleles->emplace_back(Allele::create(insertionBases, length, false));
						}
					}
					if (insertionAlleles->size() == 2) {
						std::string contig = refLoc->getContig();
						proposedEvents.emplace_back(
								VariantContextBuilder(sourceNameToAdd, contig, insertionStart, insertionStart,
								                      insertionAlleles).make());
					}
				}
				alignmentPos += elementLength;
				break;

			}
			case S: {
				alignmentPos += elementLength;
				break;
			}
			case D: {
				if (refPos > 0) {
					int deletionBasesLength;
					std::shared_ptr<uint8_t[]> deletionBases = Mutect2Utils::copyOfRange(ref, refLength, refPos - 1,
					                                                                     refPos + elementLength,
					                                                                     deletionBasesLength);
					std::shared_ptr<std::vector<std::shared_ptr<Allele>>> deletionAlleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
					int deletionStart = refLoc->getStart() + refPos - 1;
					uint8_t refByte = ref.get()[refPos - 1];
					if (BaseUtils::isRegularBase(refByte) &&
					    BaseUtils::isAllRegularBases(deletionBases, deletionBasesLength)) {
						deletionAlleles->emplace_back(Allele::create(deletionBases, deletionBasesLength, true));
						deletionAlleles->emplace_back(Allele::create(refByte, false));
						std::string contig = refLoc->getContig();
						proposedEvents.emplace_back(VariantContextBuilder(sourceNameToAdd, contig, deletionStart,
						                                                  deletionStart + elementLength,
						                                                  deletionAlleles).make());
					}
				}
				refPos += elementLength;
				break;
			}
			case M:
			case EQ:
			case X: {
				std::deque<int> mismatchOffsets;
				for (int offset = 0; offset < elementLength; offset++) {
					uint8_t refByte = BaseUtils::getUpper(ref.get()[refPos + offset]);
					uint8_t altByte = BaseUtils::getUpper(alignment.get()[alignmentPos + offset]);
					bool mismatch = refByte != altByte && BaseUtils::isRegularBase(refByte) &&
					                BaseUtils::isRegularBase(altByte);
					if (mismatch) {
						mismatchOffsets.push_back(offset);
					}
				}
				while (!mismatchOffsets.empty()) {
					int start = mismatchOffsets.front();
					mismatchOffsets.pop_front();
					int end = start;
					while (!mismatchOffsets.empty() && mismatchOffsets.front() - end <= maxMnpDistance) {
						end = mismatchOffsets.front();
						mismatchOffsets.pop_front();
					}
					int length;
					std::shared_ptr<uint8_t[]> bases = Mutect2Utils::copyOfRange(ref, refLength, refPos + start,
					                                                             refPos + end + 1, length);
					std::shared_ptr<Allele> refAllele = Allele::create(bases, length, true);
					bases = Mutect2Utils::copyOfRange(alignment, alignmentLength, alignmentPos + start,
					                                  alignmentPos + end + 1, length);
					std::shared_ptr<Allele> altAllele = Allele::create(bases, length, false);
					std::string contig = refLoc->getContig();
					std::shared_ptr<std::vector<std::shared_ptr<Allele>>> Alleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
					Alleles->emplace_back(refAllele);
					Alleles->emplace_back(altAllele);
					proposedEvents.emplace_back(
							VariantContextBuilder(sourceNameToAdd, contig, refLoc->getStart() + refPos + start,
							                      refLoc->getStart() + refPos + end, Alleles).make());
				}
				refPos += elementLength;
				alignmentPos += elementLength;
				break;
			}
			case N:
			case H:
			case P:
			default:
				throw std::invalid_argument("Unsupported cigar operator created during SW alignment");
		}
	}
	for (const std::shared_ptr<VariantContext> &proposedEvent: proposedEvents)
		addVC(proposedEvent, true);
}

void EventMap::addVC(const std::shared_ptr<VariantContext> &vc, bool merge) {
	Mutect2Utils::validateArg(vc.get(), "null is not allowed here");
	if (variantMap.find(vc->getStart()) != variantMap.end()) {
		Mutect2Utils::validateArg(merge, "Will not merge previously bound variant contexts as merge is false");
		std::shared_ptr<VariantContext> prev = variantMap.at(vc->getStart());
		variantMap[vc->getStart()] = makeBlock(prev, vc);
	} else
		variantMap.insert(std::make_pair(vc->getStart(), vc));
}

std::shared_ptr<VariantContext>
EventMap::makeBlock(std::shared_ptr<VariantContext> vc1, std::shared_ptr<VariantContext> vc2) {
	Mutect2Utils::validateArg(vc1->getStart() == vc2->getStart(), "vc1 and 2 must have the same start");
	Mutect2Utils::validateArg(vc1->isBiallelic(), "vc1 must be biallelic");
	if (!vc1->isSNP()) {
		Mutect2Utils::validateArg((vc1->isSimpleDeletion() && vc2->isSimpleInsertion()) ||
		                          (vc1->isSimpleInsertion() && vc2->isSimpleDeletion()),
		                          "Can only merge single insertion with deletion (or vice versa)");
	} else {
		Mutect2Utils::validateArg(!vc2->isSNP(), "there's been some terrible bug in the cigar");
	}

	std::shared_ptr<Allele> new_ref;
	std::shared_ptr<Allele> new_alt;
	VariantContextBuilder b(vc1);
	if (vc1->isSNP()) {
		if (*(vc1->getReference()) == (*(vc2->getReference()))) {
			new_ref = vc1->getReference();
			std::shared_ptr<uint8_t[]> vc1Bases = vc1->getAlternateAllele(0)->getBases();
			int vc1Length = vc1->getAlternateAllele(0)->getBasesLength();
			std::shared_ptr<uint8_t[]> vc2Bases = vc2->getAlternateAllele(0)->getBases();
			int vc2Length = vc2->getAlternateAllele(0)->getBasesLength();
			int tmpLength = vc1Length + vc2Length - 1;
			std::shared_ptr<uint8_t[]> tmp(new uint8_t[tmpLength]);
			memcpy(tmp.get(), vc1Bases.get(), vc1Length);
			memcpy(tmp.get() + vc1Length, vc2Bases.get() + 1, vc2Length - 1);
			new_alt = Allele::create(tmp, tmpLength, false);
		} else {
			new_ref = vc2->getReference();
			new_alt = vc1->getAlternateAllele(0);
			b.setStop(vc2->getEnd());
		}
	} else {
		std::shared_ptr<VariantContext> insertion = vc1->isSimpleInsertion() ? vc1 : vc2;
		std::shared_ptr<VariantContext> deletion = vc1->isSimpleInsertion() ? vc2 : vc1;
		new_ref = deletion->getReference();
		new_alt = insertion->getAlternateAllele(0);
		b.setStop(deletion->getEnd());
	}
	std::shared_ptr<std::vector<std::shared_ptr<Allele>>> Alleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
	Alleles->emplace_back(new_ref);
	Alleles->emplace_back(new_alt);
	b.setAlleles(Alleles);
	return b.make();
}

bool EventMap::empty() {
	return variantMap.empty();
}

std::set<int> EventMap::buildEventMapsForHaplotypes(std::vector<std::shared_ptr<Haplotype>> &haplotypes,
                                                    const std::shared_ptr<uint8_t[]> &ref, int refLength,
                                                    const std::shared_ptr<Locatable> &refLoc, int maxMnpDistance) {
	Mutect2Utils::validateArg(maxMnpDistance >= 0, "maxMnpDistance may not be negative.");
	std::set<int> startPosKeySet;
	int hapNumber = 0;
	for (auto &h: haplotypes) {
		auto *newEventMap = new EventMap(h, ref, refLength, refLoc, "HC" + std::to_string(hapNumber), maxMnpDistance);
        // deal with memory leak
        if (h->getEventMap() != nullptr)
            delete h->getEventMap();
		h->setEventMap(newEventMap);
		/*if (h->getIsNonReference()) {
			std::cout.precision(10);
			std::cout.flags(std::ostream::fixed);
			std::cout << h->getScore() << "\t" << h->getBasesLength() << std::endl;
			std::cout << h->getBaseString() << std::endl;
			std::cout << "variantMap size " << newEventMap->variantMap.size() << std::endl;
			for (const auto &ve: newEventMap->variantMap) {
				std::cout << ve.second->getStart() + 1 << " " << ve.second->getEnd() + 1 << std::endl;
			}
		}*/
		for (int i: newEventMap->getStartPositions()) {
			startPosKeySet.insert(i);
		}
	}
	return startPosKeySet;
}

std::set<int> EventMap::getStartPositions() {
	std::set<int> ret;
	for (std::pair<int, std::shared_ptr<VariantContext>> pair_vc: variantMap) {
		ret.insert(pair_vc.first);
	}
	return {ret.begin(), ret.end()};
}

std::set<std::shared_ptr<VariantContext>, VariantContextComparator>
EventMap::getAllVariantContexts(std::vector<std::shared_ptr<Haplotype>> &haplotypes) {
	std::set<std::shared_ptr<VariantContext>, VariantContextComparator> ret;
	for (const std::shared_ptr<Haplotype> &h: haplotypes) {
		for (const std::shared_ptr<VariantContext> &vc: h->getEventMap()->getVariantContexts()) {
			ret.insert(vc);
		}
	}
	return ret;
}

std::vector<std::shared_ptr<VariantContext>> EventMap::getVariantContexts() {
	std::vector<std::shared_ptr<VariantContext>> ret;
	for (std::pair<int, std::shared_ptr<VariantContext>> pair_vc: variantMap) {
		ret.emplace_back(pair_vc.second);
	}
	return ret;
}

std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> EventMap::getOverlappingEvents(int loc)
{
    std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> events = std::make_shared<std::vector<std::shared_ptr<VariantContext>>>();
    for(auto& iter : variantMap)
    {
        if(iter.first <= loc && iter.second->getEnd() >= loc)
        {
            events->emplace_back(iter.second);
        }
    }
    return events;
}
