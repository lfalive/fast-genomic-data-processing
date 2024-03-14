//
// Created by 梦想家xixi on 2022/1/15.
//

#include "ReferenceConfidenceModel.h"

#include <utility>

std::shared_ptr<Haplotype>
ReferenceConfidenceModel::createReferenceHaplotype(const std::shared_ptr<AssemblyRegion> &activeRegion,
                                                   const std::shared_ptr<uint8_t[]>& refBases, int &length,
                                                   const std::shared_ptr<SimpleInterval> &paddedReferenceLoc) {
	int alignmentStart = activeRegion->getExtendedSpan()->getStart() - paddedReferenceLoc->getStart();
	if (alignmentStart < 0)
		throw std::invalid_argument("Bad alignment start in createReferenceHaplotype");
	std::shared_ptr<Haplotype> refHaplotype = std::make_shared<Haplotype>(refBases, length, true);
	refHaplotype->setAlignmentStartHapwrtRef(alignmentStart);
	std::shared_ptr<Cigar> c(new Cigar());
	c->add(CigarElement(refHaplotype->getBasesLength(), M));
	refHaplotype->setCigar(c);
	return refHaplotype;
}
