//
// Created by 梦想家xixi on 2021/12/1.
//

#include "AssemblyResultSet.h"

#include <utility>
#include <algorithm>

bool AssemblyResultSet::add(const std::shared_ptr<Haplotype> &h, const std::shared_ptr<AssemblyResult> &ar) {
	Mutect2Utils::validateArg(h.get(), "input haplotype cannot be null");
	Mutect2Utils::validateArg(ar.get(), "input assembly-result cannot be null");
	Mutect2Utils::validateArg(h->getGenomeLocation().get(), "the haplotype provided must have a genomic location");

	bool assemblyResultAdditionReturn = add(ar);

	if (haplotypes.find(h) != haplotypes.end()) {
		if (assemblyResultByHaplotype.find(h) == assemblyResultByHaplotype.end()) {
			assemblyResultByHaplotype.insert(std::make_pair(h, ar));
			return true;
		}
		if (assemblyResultByHaplotype.at(h) != ar)
			throw std::invalid_argument("there is already a different assembly result for the input haplotype");
		return assemblyResultAdditionReturn;
	}
	haplotypes.insert(h);
	assemblyResultByHaplotype.insert(std::make_pair(h, ar));
	updateReferenceHaplotype(h);
	if (h->getIsNonReference()) {
		variationPresent = true;
	}
	return true;
}

bool AssemblyResultSet::add(const std::shared_ptr<AssemblyResult> &ar) {
	Mutect2Utils::validateArg(ar.get(), "input assembly-result cannot be null");
	int kmerSize = ar->getKmerSize();
	if (assemblyResultByKmerSize.find(kmerSize) != assemblyResultByKmerSize.end()) {
		if (assemblyResultByKmerSize.at(kmerSize) != ar)
			throw std::invalid_argument("a different assembly result with the same kmerSize was already added");
		return false;
	}
	assemblyResultByKmerSize.insert(std::make_pair(kmerSize, ar));
	kmerSizes.insert(kmerSize);
	return true;
}

void AssemblyResultSet::updateReferenceHaplotype(const std::shared_ptr<Haplotype> &newHaplotype) {
	if (!newHaplotype->getIsReference())
		return;
	if (refHaplotype == nullptr) {
		refHaplotype = newHaplotype;
	} else
		throw std::invalid_argument("the assembly-result-set already have a reference haplotype that is different");
}

void AssemblyResultSet::setRegionForGenotyping(std::shared_ptr<AssemblyRegion> regionForGenotyping) {
	this->regionForGenotyping = regionForGenotyping;
}

void AssemblyResultSet::setFullReferenceWithPadding(std::shared_ptr<uint8_t[]> fullReferenceWithPadding, int length) {
	this->fullReferenceWithPadding = fullReferenceWithPadding;
	this->fullReferenceWithPaddingLength = length;
}

void AssemblyResultSet::setPaddedReferenceLoc(const std::shared_ptr<SimpleInterval> &paddedReferenceLoc) {
	this->paddedReferenceLoc = paddedReferenceLoc;
}

bool AssemblyResultSet::add(const std::shared_ptr<Haplotype> &h) {
	Mutect2Utils::validateArg(h.get(), "input haplotype can not be null");
	Mutect2Utils::validateArg(h->getGenomeLocation().get(), "haplotype genomeLocation cannot be null");
	if (haplotypes.find(h) != haplotypes.end())
		return false;
	haplotypes.insert(h);
	updateReferenceHaplotype(h);
	return true;
}

std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &
AssemblyResultSet::getVariationEvents(int maxMnpDistance) {
	if (maxMnpDistance < 0)
		throw std::invalid_argument("maxMnpDistance may not be negative.");
	bool sameMnpDistance = lastMaxMnpDistanceUsed != -1 && maxMnpDistance == lastMaxMnpDistanceUsed;
	lastMaxMnpDistanceUsed = maxMnpDistance;
	bool flag = false;
	for (const auto &haplotype: haplotypes) {
		if (haplotype->getIsNonReference() &&
		    (haplotype->getEventMap() == nullptr || haplotype->getEventMap()->empty())) {
			flag = true;
			break;
		}
	}
	if (variationEvents.empty() || !sameMnpDistance || flag) {
		regenerateVariationEvents(maxMnpDistance);
	}
	return variationEvents;
}

void AssemblyResultSet::regenerateVariationEvents(int maxMnpDistance) {
	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> haplotypeList =
#ifdef SORT_MODE
			getSortedHaplotypeList();
#else
			getHaplotypeList();
#endif
	EventMap::buildEventMapsForHaplotypes(*haplotypeList, fullReferenceWithPadding, fullReferenceWithPaddingLength,
	                                      paddedReferenceLoc, maxMnpDistance);
	variationEvents = EventMap::getAllVariantContexts(*haplotypeList);
	lastMaxMnpDistanceUsed = maxMnpDistance;
	for (const std::shared_ptr<Haplotype> &haplotype: *haplotypeList) {
		if (haplotype->getIsNonReference()) {
			variationPresent = true;
			break;
		}
	}
}

std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>>
AssemblyResultSet::getSortedHaplotypeList() {
	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>>
			res = std::make_shared<std::vector<std::shared_ptr<Haplotype>>>();
	res->reserve(haplotypes.size());
	for (const auto &haplotype: haplotypes) {
		res->emplace_back(haplotype);
	}
	std::sort(res->begin(), res->end(),
	          [](const std::shared_ptr<Haplotype> &h1, const std::shared_ptr<Haplotype> &h2) -> bool {
		          double score1 = h1->getScore(), score2 = h2->getScore();
		          if (score1 - score2 > 0.0000000001)
			          return true;
		          if (score2 - score1 > 0.0000000001)
			          return false;
		          int len1 = h1->getBasesLength(), len2 = h2->getBasesLength();
		          if (len1 != len2)
			          return len1 > len2;
		          std::shared_ptr<uint8_t[]> base1 = h1->getBases(), base2 = h2->getBases();
		          for (int i = 0; i < len1; ++i) {
			          if (base1[i] == base2[i]) continue;
			          return base1[i] > base2[i];
		          }
		          return false;
	          });
	return res;
}

std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>>
AssemblyResultSet::getHaplotypeList() {
	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>>
			res = std::make_shared<std::vector<std::shared_ptr<Haplotype>>>();
	res->reserve(haplotypes.size());
	for (const auto &haplotype: haplotypes) {
		res->emplace_back(haplotype);
	}
	return res;
}

std::vector<std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>>>
AssemblyResultSet::calculateOriginalByTrimmedHaplotypes(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion) {
	std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype >>
			originalByTrimmedHaplotypes =
#ifdef SORT_MODE
			trimDownHaplotypes(trimmedAssemblyRegion, getSortedHaplotypeList());
#else
			trimDownHaplotypes(trimmedAssemblyRegion, getHaplotypeList());
#endif
	std::vector<std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>>> sortedOriginalByTrimmedHaplotypes;
	sortedOriginalByTrimmedHaplotypes.reserve(originalByTrimmedHaplotypes->size());
	for (const auto &element: *originalByTrimmedHaplotypes) {
		sortedOriginalByTrimmedHaplotypes.emplace_back(std::make_pair(element.first, element.second));
	}
	std::sort(sortedOriginalByTrimmedHaplotypes.begin(), sortedOriginalByTrimmedHaplotypes.end(),
	          [](const std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>> &left,
	             const std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>> &right) {
		          std::shared_ptr<Haplotype> ha = left.first, hb = right.first;
		          int len1 = ha->getBasesLength(), len2 = hb->getBasesLength();
		          if (len1 != len2)
			          return len1 < len2;
		          uint8_t *seq1 = ha->getBases().get(), *seq2 = hb->getBases().get();
		          for (int k = 0; k < len1; ++k) {
			          if (seq1[k] == seq2[k]) continue;
			          return seq1[k] < seq2[k];
		          }
		          return false;
	          });
	return sortedOriginalByTrimmedHaplotypes;
}

std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
AssemblyResultSet::trimDownHaplotypes(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion,
                                      const std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> &haplotypeList) {
	std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>
			originalByTrimmedHaplotypes = std::make_shared<phmap::flat_hash_map<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype>>();
	for (const auto &h: *haplotypeList) {
		std::shared_ptr<Haplotype> trimmed = h->trim(trimmedAssemblyRegion->getExtendedSpan());
		if (trimmed != nullptr) {
			auto iter = originalByTrimmedHaplotypes->find(trimmed);
			if (iter != originalByTrimmedHaplotypes->end()) {
				if (trimmed->getIsReference()) {
					originalByTrimmedHaplotypes->erase(iter);
					originalByTrimmedHaplotypes->insert({trimmed, h});
				}
			} else {
				originalByTrimmedHaplotypes->insert({trimmed, h});
			}
		} else if (h->getIsReference())
			throw std::invalid_argument("trimming eliminates the reference haplotype");
	}
	return originalByTrimmedHaplotypes;
}

std::shared_ptr<AssemblyResultSet>
AssemblyResultSet::trimTo(const std::shared_ptr<AssemblyRegion> &trimmedAssemblyRegion) {
	std::vector<std::pair<std::shared_ptr<Haplotype>, std::shared_ptr<Haplotype>>>
			originalByTrimmedHaplotypes = calculateOriginalByTrimmedHaplotypes(trimmedAssemblyRegion);
	if (refHaplotype == nullptr)
		throw std::invalid_argument("refHaplotype is null");

	std::shared_ptr<AssemblyResultSet> result = std::make_shared<AssemblyResultSet>();
	for (const auto &element: originalByTrimmedHaplotypes) {
		const std::shared_ptr<Haplotype> &trimmed = element.first;
		const std::shared_ptr<Haplotype> &original = element.second;
		if (original == nullptr)
			throw std::invalid_argument("all trimmed haplotypes must have an original one");

		if (assemblyResultByHaplotype.find(original) == assemblyResultByHaplotype.end())
			result->add(trimmed);
		else {
			std::shared_ptr<AssemblyResult> &as = assemblyResultByHaplotype.at(original);
			result->add(trimmed, as);
		}
	}

	result->setRegionForGenotyping(trimmedAssemblyRegion);
	result->setFullReferenceWithPadding(fullReferenceWithPadding, fullReferenceWithPaddingLength);
	result->setPaddedReferenceLoc(paddedReferenceLoc);
	result->variationPresent = false;
	for (const auto &haplotype: haplotypes) {
		if (haplotype->getIsNonReference()) {
			result->variationPresent = true;
			break;
		}
	}
	result->wasTrimmed = true;
	return result;
}

bool AssemblyResultSet::isisVariationPresent() {
	return variationPresent && haplotypes.size() > 1;
}

std::shared_ptr<AssemblyRegion> AssemblyResultSet::getRegionForGenotyping() {
	return regionForGenotyping;
}

std::shared_ptr<Haplotype> &AssemblyResultSet::getReferenceHaplotype() {
	return refHaplotype;
}

std::shared_ptr<SimpleInterval> &AssemblyResultSet::getPaddedReferenceLoc() {
	return paddedReferenceLoc;
}

std::shared_ptr<uint8_t[]> &AssemblyResultSet::getFullReferenceWithPadding() {
	return fullReferenceWithPadding;
}

int AssemblyResultSet::getFullReferenceWithPaddingLength() {
	return fullReferenceWithPaddingLength;
}

void AssemblyResultSet::deleteEventMap() {
	auto haplotypesToReleased = *this->getHaplotypeList();
	for (auto &item: haplotypesToReleased) {
		if (item->getEventMap() != nullptr) {
			delete item->getEventMap();
			item->setEventMap(nullptr);
		}
	}
}

void AssemblyResultSet::printSortedHaplotypes() {
	std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> ret = getSortedHaplotypeList();
	std::cout << "Haplotypes\t" << ret->size() << std::endl;
	for (const auto &h: *ret) {
		std::string baseStr = h->getBaseString();
		if (h->getIsNonReference()) {
			std::cout.precision(10);
			std::cout.flags(std::ostream::fixed);
			std::cout << h->getScore() << "\t";
		} else {
			std::cout << "ref\t";
		}
		std::cout << baseStr.length() << " ";
		for (const auto &ce: h->getCigar()->getCigarElements()) {
			std::cout << ce.getLength() << CigarOperatorUtils::enumToCharacter(ce.getOperator());
		}
		std::cout << std::endl << baseStr << std::endl;;
	}
}
