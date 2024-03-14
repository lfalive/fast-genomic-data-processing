//
// Created by 梦想家xixi on 2021/11/10.
//

#include "AlignmentUtils.h"
#include <memory>
#include <utility>
#include <cstring>
#include <cassert>
#include "CigarUtils.h"
#include "Mutect2Utils.h"

std::vector<CigarPairTransform> AlignmentUtils::cigarPairTransformers = std::vector<CigarPairTransform>{
    //
    // op12 is a match
    //
    // 3: xxx B yyy
    // ^^^^^^^^^^^^
    // 2: xxx M yyy
    // 1: xxx M yyy
    CigarPairTransform(CigarOperator::M, CigarOperator::M, CigarOperator::M, 1, 1),
    // 3: xxx I yyy
    // ^^^^^^^^^^^^
    // 2: xxx I yyy
    // 1: xxx M yyy
    CigarPairTransform(CigarOperator::M, CigarOperator::I, CigarOperator::I, 1, 1),
    // 3: xxx D yyy
    // ^^^^^^^^^^^^
    // 2: xxx D yyy
    // 1: xxx M yyy
    CigarPairTransform(CigarOperator::M, CigarOperator::D, CigarOperator::D, 0, 1),

    //
    // op12 is a deletion
    //
    // 3: xxx D M yyy
    // ^^^^^^^^^^^^
    // 2: xxx M yyy
    // 1: xxx D yyy
    CigarPairTransform(CigarOperator::D, CigarOperator::M, CigarOperator::D, 1, 1),
    // 3: xxx D2 D1 yyy
    // ^^^^^^^^^^^^
    // 2: xxx D2 yyy
    // 1: xxx D1 yyy
    CigarPairTransform(CigarOperator::D, CigarOperator::D, CigarOperator::D, 0, 1),
    // 3: xxx X yyy => no-op, we skip emitting anything here
    // ^^^^^^^^^^^^
    // 2: xxx I yyy
    // 1: xxx D yyy
    CigarPairTransform(CigarOperator::D, CigarOperator::I, CigarOperator::_NULL, 1, 1),

    //
    // op12 is a insertion
    //
    // 3: xxx I M yyy
    // ^^^^^^^^^^^^
    // 2: xxx M yyy
    // 1: xxx I yyy
    CigarPairTransform(CigarOperator::I, CigarOperator::M, CigarOperator::I, 1, 0),
    // 3: xxx I D yyy
    // ^^^^^^^^^^^^
    // 2: xxx D yyy
    // 1: xxx I yyy
     CigarPairTransform(CigarOperator::I, CigarOperator::D, CigarOperator::I, 1, 0),
    // 3: xxx I1 I2 yyy
    // ^^^^^^^^^^^^
    // 2: xxx I2 yyy
    // 1: xxx I1 yyy
     CigarPairTransform(CigarOperator::I, CigarOperator::I, CigarOperator::I, 1, 0)
};

std::set<CigarOperator> AlignmentUtils::ALIGNED_TO_GENOME_PLUS_SOFTCLIPS = {CigarOperator::M, CigarOperator::EQ, CigarOperator::X, CigarOperator::S};

std::shared_ptr<Cigar> AlignmentUtils::consolidateCigar(std::shared_ptr<Cigar> c) {
    if(c == nullptr) {throw std::invalid_argument("Cigar cannot be null");}

    if(!needsConsolidation(c))
        return c;

    std::shared_ptr<Cigar> returnCigar(new Cigar());
    int sumLength = 0;
    CigarElement* lastElement = nullptr;

    for(CigarElement & cur : c->getCigarElements()) {
        if( cur.getLength() == 0)
            continue;

        if(lastElement != nullptr && lastElement->getOperator() != cur.getOperator()) {
            returnCigar->add(CigarElement(sumLength, lastElement->getOperator()));
            sumLength = 0;
        }
        sumLength += cur.getLength();
        lastElement = &cur;
    }
    if(sumLength > 0) {
        returnCigar->add(CigarElement(sumLength, lastElement->getOperator()));
    }
    return returnCigar;
}

bool AlignmentUtils::needsConsolidation(const std::shared_ptr<Cigar>& c) {
    if(c->numCigarElements() <= 1)
        return false;
    CigarOperator lastOp;
    for(CigarElement cur : c->getCigarElements()) {
        if(cur.getLength() == 0 || lastOp == cur.getOperator())
            return true;
        lastOp = cur.getOperator();
    }
    return false;
}

std::pair<int, std::shared_ptr<uint8_t[]>> AlignmentUtils::getBasesCoveringRefInterval(int refStart, int refEnd, std::shared_ptr<uint8_t[]>bases, int length, int basesStartOnRef,
                                                     const std::shared_ptr<Cigar>& basesToRefCigar) {
    if(refStart < 0 || refEnd < refStart) {
        throw std::invalid_argument("Bad start and/or stop");
    }
    Mutect2Utils::validateArg(basesStartOnRef >= 0, "BasesStartOnRef must be >= 0");
    Mutect2Utils::validateArg(bases != nullptr, "Null is not allowed there");
    Mutect2Utils::validateArg(basesToRefCigar != nullptr, "Null is not allowed there");
    Mutect2Utils::validateArg(length == basesToRefCigar->getReadLength(), "Mismatch in length between reference bases and cigar length.");

    int refPos = basesStartOnRef;
    int basesPos = 0;
    int basesStart = -1;
    int basesStop = -1;
    bool done = false;

    for(int iii = 0; ! done && iii <basesToRefCigar->numCigarElements(); iii++) {
        CigarElement ce = basesToRefCigar->getCigarElement(iii);
        switch (ce.getOperator()) {
            case I:
                basesPos += ce.getLength();
                break;
            case M:
            case X:
            case EQ:
                for (int i = 0; i < ce.getLength(); i++) {
                    if(refPos == refStart)
                        basesStart = basesPos;
                    if(refPos == refEnd) {
                        basesStop = basesPos;
                        done = true;
                        break;
                    }
                    refPos++;
                    basesPos++;
                }
                break;
            case D:
                for(int i = 0; i < ce.getLength(); i++) {
                    if(refPos == refEnd || refPos == refStart) {
                        return {0, nullptr};
                    }
                    refPos++;
                }
                break;
            default:
                throw std::invalid_argument("Unsupported operator");
        }
    }

    if(basesStart == -1 || basesStop == -1)
        throw std::invalid_argument("Never found start or stop");

    int newLength = basesStop - basesStart + 1;
    std::shared_ptr<uint8_t[]> ret(new uint8_t[newLength]);
    memcpy(ret.get(), bases.get() + basesStart, newLength);
    return {newLength, ret};
}

std::shared_ptr<Cigar> AlignmentUtils::trimCigarByReference(const std::shared_ptr<Cigar>&cigar, const int start, const int end) {
    Mutect2Utils::validateArg(start >= 0, "Start must be >= 0");
    Mutect2Utils::validateArg(end >= start, "End is < start");
    Mutect2Utils::validateArg(end <= cigar->getReferenceLength(), "End is beyond the cigar's reference length");

    std::shared_ptr<Cigar> result = trimCigar(cigar, start, end, true);
    Mutect2Utils::validateArg(result->getReferenceLength() == end - start + 1, "trimCigarByReference failure");
    return result;
}

std::shared_ptr<Cigar> AlignmentUtils::trimCigar(const std::shared_ptr<Cigar>& cigar, const int start, const int end, const bool byReference) {
    std::vector<CigarElement> newElements;

    int pos = 0;
    for (CigarElement elt : cigar->getCigarElements()) {
        if ( pos > end && (byReference || elt.getOperator() != D) ) break;

        switch (elt.getOperator()) {
            case D:
                if(!byReference) {
                    if(pos >= start)
                        newElements.emplace_back(elt);
                    break;
                }
            case EQ:
            case M:
            case X:
                pos = addCigarElements(newElements, pos, start, end, elt);
                break;
            case S:
            case I:
                if(byReference) {
                    if(pos >= start)
                        newElements.emplace_back(elt);
                } else {
                    pos = addCigarElements(newElements, pos, start, end, elt);
                }
                break;
            default:
                throw std::invalid_argument("Cannot handle");
        }
    }
    std::shared_ptr<Cigar> tmp(new Cigar(newElements));
    std::shared_ptr<Cigar> ret = consolidateCigar(tmp);
    return ret;
}

int AlignmentUtils::addCigarElements(std::vector<CigarElement> & dest, int pos, int start, int end, CigarElement elt) {
    int length = std::min(pos + elt.getLength() - 1, end) - std::max(pos, start) + 1;
    if(length > 0)
        dest.emplace_back(CigarElement(length, elt.getOperator()));
    return pos + elt.getLength();
}

bool AlignmentUtils::startsOrEndsWithInsertionOrDeletion(const std::shared_ptr<Cigar>& cigar) {
    Mutect2Utils::validateArg(cigar != nullptr, "Null is not allowed");
    if(cigar->isEmpty())
        return false;
    CigarOperator first = cigar->getCigarElement(0).getOperator();
    CigarOperator last = cigar->getCigarElement(cigar->numCigarElements()-1).getOperator();
    return first == D || first == I || last == D || last == I;
}

std::shared_ptr<Cigar> AlignmentUtils::removeTrailingDeletions(std::shared_ptr<Cigar> c) {
    std::vector<CigarElement> elements = c->getCigarElements();
    if(elements.at(elements.size()-1).getOperator() != D)
        return c;
    std::vector<CigarElement> newElements(elements.begin(), elements.end()-1);
    return std::make_shared<Cigar>(newElements);
}

std::shared_ptr<Cigar> AlignmentUtils::trimCigarByBases(const std::shared_ptr<Cigar>&cigar, int start, int end) {
    Mutect2Utils::validateArg(start >= 0, "tart must be >= 0 but got");
    Mutect2Utils::validateArg(end >= start, "End is < start ");
    Mutect2Utils::validateArg(end <= (cigar->getReadLength()), "End is beyond the cigar's read length");

    std::shared_ptr<Cigar> result = trimCigar(cigar, start, end, false);
    int expectedSize = end - start + 1;
    Mutect2Utils::validateArg((result->getReadLength()) == expectedSize, "trimCigarByBases failure");
    return result;
}

std::shared_ptr<Cigar>
AlignmentUtils::leftAlignIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength,
                               std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex,
                               int leftmostAllowedAlignment, bool doNotThrowExceptionForMultipleIndels) {
    ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

    int numIndels = countIndelElements(cigar);
    if ( numIndels == 0 )
        return cigar;
    if ( numIndels == 1 )
        return leftAlignSingleIndel(cigar, refSeq, refLength, readSeq, readLength, refIndex, readIndex, leftmostAllowedAlignment, true);

    // if we got here then there is more than 1 indel in the alignment
    if ( doNotThrowExceptionForMultipleIndels )
        return cigar;

    throw "attempting to left align a CIGAR that has more than 1 indel in its alignment but this functionality has not been implemented yet";
}

std::shared_ptr<Cigar>
AlignmentUtils::leftAlignIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength,
                               std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex, bool doNotThrowExceptionForMultipleIndels){
    return leftAlignIndel(cigar, refSeq, refLength, readSeq, readLength, refIndex, readIndex, 0, doNotThrowExceptionForMultipleIndels);
}

int AlignmentUtils::countIndelElements(std::shared_ptr<Cigar> &cigar)
{
    int indelCount = 0;
    for ( CigarElement& ce : cigar->getCigarElements() ) {
        if ( ce.getOperator() == CigarOperator::D || ce.getOperator() == CigarOperator::I )
            indelCount++;
    }
    return indelCount;
}

std::shared_ptr<Cigar>
AlignmentUtils::leftAlignSingleIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]>refSeq, int refLength, std::shared_ptr<uint8_t[]>readSeq, int readLength,
                                     int refIndex, int readIndex, bool cleanupCigar) {
    return leftAlignSingleIndel(std::move(cigar), refSeq, refLength, readSeq, readLength, refIndex, readIndex, 0, false);
}

std::shared_ptr<Cigar>
AlignmentUtils::leftAlignSingleIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]>refSeq, int refLength, std::shared_ptr<uint8_t[]>readSeq, int readLength,
                                     int refIndex, int readIndex, int leftmostAllowedAlignment, bool cleanupCigar1) {
    ensureLeftAlignmentHasGoodArguments(cigar, refSeq, readSeq, refIndex, readIndex);

    int indexOfIndel = -1;
    for(int i = 0; i < cigar->numCigarElements(); i++) {
        CigarElement ce = cigar->getCigarElement(i);
        if(ce.getOperator() == D || ce.getOperator() == I) {
            if(indexOfIndel != -1)
                throw std::invalid_argument("attempting to left align a CIGAR that has more than 1 indel in its alignment");
            indexOfIndel = i;
        }
    }

    if(indexOfIndel == -1)
        throw std::invalid_argument("attempting to left align a CIGAR that has no indels in its alignment");
    if(indexOfIndel == 0)
        return cigar;

    int indelLength = cigar->getCigarElement(indexOfIndel).getLength();
    int altStringLength = 0;
	/*std::cout<<std::string ((char*)refSeq.get(), refLength)<<std::endl;
	std::cout<<std::string ((char*)readSeq.get(), readLength)<<std::endl;
	for (const auto &item: cigar->getCigarElements()){
		std::cout<<item.getLength()<<CigarOperatorUtils::enumToCharacter(item.getOperator());
	}
	std::cout<<" refLength:"<<refLength<<" readLength:"<<readLength<<" refIndex:"<<refIndex<<" readIndex:"<<readIndex<<std::endl;*/
	std::shared_ptr<uint8_t[]> altString = createIndelString(cigar, indexOfIndel, refSeq, refLength, readSeq, readLength, refIndex, readIndex, altStringLength);
    if(altString == nullptr)
        return cigar;
    std::shared_ptr<Cigar> newCigar(new Cigar(*cigar));
    for(int i = 0; i < indelLength; i++) {
        std::shared_ptr<Cigar> tmp = moveCigarLeft(newCigar, indexOfIndel);
        newCigar = tmp;
        if(isIndelAlignedTooFarLeft(newCigar,leftmostAllowedAlignment)){
            break;
        }
        int newAltStringLength = 0;
        std::shared_ptr<uint8_t[]> newAltString = createIndelString(newCigar, indexOfIndel, refSeq, refLength, readSeq, readLength, refIndex, readIndex, newAltStringLength);
        bool reachedEndOfRead = cigarHasZeroSizeElement(newCigar);
        bool isEqual = true;
        if(altStringLength != newAltStringLength)
            isEqual = false;
        else {
            uint8_t * alt = altString.get();
            uint8_t * newAlt = newAltString.get();
            for(int j = 0; j < altStringLength; j++) {
                if(alt[j] != newAlt[j]){
                    isEqual = false;
                    break;
                }
            }
        }
        if(isEqual) {
            cigar = newCigar;
            i = -1;
            if(reachedEndOfRead)
                cigar = cleanupCigar1 ? cleanUpCigar(cigar) : cigar;
        }
        if(reachedEndOfRead)
            break;
    }
    return cigar;
}

void AlignmentUtils::ensureLeftAlignmentHasGoodArguments(const std::shared_ptr<Cigar>&cigar, std::shared_ptr<uint8_t[]>refSeq, std::shared_ptr<uint8_t[]>readSeq,
                                                          int refIndex, int readIndex) {
    Mutect2Utils::validateArg(cigar.get(), "ciagr");
    Mutect2Utils::validateArg(refSeq.get(), "refSeq");
    Mutect2Utils::validateArg(readSeq.get(), "readSeq");
    Mutect2Utils::validateArg(refIndex >= 0, "attempting to left align with a reference index less than 0");
    Mutect2Utils::validateArg(readIndex >= 0, "attempting to left align with a read index less than 0");
}

std::shared_ptr<uint8_t[]>
AlignmentUtils::createIndelString(const std::shared_ptr<Cigar>& cigar, int indexOfIndel, std::shared_ptr<uint8_t[]>refSeq, int refLength, std::shared_ptr<uint8_t[]>readSeq,
                                  int readLength, int refIndex, int readIndex, int & newLength) {
    CigarElement indel = cigar->getCigarElement(indexOfIndel);
    int indelLength = indel.getLength();

    int totalRefBases = 0;
    for(int i = 0; i < indexOfIndel; i++) {
        CigarElement ce = cigar->getCigarElement(i);
        int length = ce.getLength();

        switch (ce.getOperator()) {
            case M:
            case EQ:
            case X:
                readIndex += length;
                refIndex += length;
                totalRefBases += length;
                break;
            case S:
                readIndex += length;
                break;
            case N:
                refIndex += length;
                totalRefBases += length;
                break;
            default:
                break;
        }
    }

    if(totalRefBases + indelLength > refLength)
        indelLength -= (totalRefBases + indelLength - refLength);

    int altLength = refLength + (indelLength * (indel.getOperator() == D ? -1 : 1));
    std::shared_ptr<uint8_t[]> alt{new uint8_t[altLength]};

    if(refIndex > altLength || refIndex > refLength)
        return nullptr;
    memcpy(alt.get(), refSeq.get(), refIndex);
    int currentPos = refIndex;

    if(indel.getOperator() == D) {
        refIndex += indelLength;
    } else {
        memcpy(alt.get()+currentPos, readSeq.get()+readIndex, indelLength);
        currentPos += indelLength;
    }

    if(refLength - refIndex > altLength - currentPos)
        return nullptr;
    memcpy(alt.get()+currentPos, refSeq.get()+refIndex, refLength-refIndex);
    newLength = altLength;
    return alt;
}

std::shared_ptr<Cigar> AlignmentUtils::moveCigarLeft(const std::shared_ptr<Cigar>&cigar, int indexOfIndel) {
    std::vector<CigarElement> elements;
    elements.reserve(indexOfIndel - 1);
    for(int i = 0; i < indexOfIndel - 1; i++)
        elements.emplace_back(cigar->getCigarElement(i));

    CigarElement ce = cigar->getCigarElement(indexOfIndel - 1);
    elements.emplace_back(CigarElement(std::max(ce.getLength()-1, 0), ce.getOperator()));
    elements.emplace_back(cigar->getCigarElement(indexOfIndel));
    if(indexOfIndel + 1 < cigar->numCigarElements()) {
        ce = cigar->getCigarElement(indexOfIndel + 1);
        elements.emplace_back(CigarElement(ce.getLength() + 1, ce.getOperator()));
    } else {
        elements.emplace_back(CigarElement(1, M));
    }

    for (int i = indexOfIndel + 2; i < cigar->numCigarElements(); i++)
        elements.emplace_back(cigar->getCigarElement(i));

    return std::make_shared<Cigar>(elements);
}

bool AlignmentUtils::isIndelAlignedTooFarLeft(const std::shared_ptr<Cigar>&cigar, int leftmostAllowedAlignment) {
    int location = 0;
    for(CigarElement element : cigar->getCigarElements()) {
        if(element.getOperator() == D || element.getOperator() == I) {
            return location<leftmostAllowedAlignment;
        }
        if(CigarOperatorUtils::getConsumesReferenceBases(element.getOperator())) {
            location += element.getLength();
        }
    }
    return false;
}

bool AlignmentUtils::cigarHasZeroSizeElement(const std::shared_ptr<Cigar>&c) {
    auto & cigarElements =  c->getCigarElements();
    if(std::any_of(cigarElements.begin(), cigarElements.end(), [](CigarElement& ce){return ce.getLength() == 0;}))
        return true;

    return false;

    for(CigarElement ce : c->getCigarElements()) {
        if(ce.getLength() == 0)
            return true;
    }
    return false;
}

std::shared_ptr<Cigar> AlignmentUtils::cleanUpCigar(const std::shared_ptr<Cigar>& c) {
    std::vector<CigarElement> elements;
    for(CigarElement ce : c->getCigarElements()) {
        if(ce.getLength() != 0 && (! elements.empty() || ce.getOperator() != D)) {
            elements.emplace_back(ce);
        }
    }

    return std::make_shared<Cigar>(elements);
}

std::shared_ptr<SAMRecord>
AlignmentUtils::createReadAlignedToRef(const std::shared_ptr<SAMRecord>& originalRead, const std::shared_ptr<Haplotype>& haplotype,
                                       const std::shared_ptr<Haplotype>& refHaplotype, int referenceStart, bool isInformative,
                                       SmithWatermanAligner *aligner) {
    assert(originalRead);
    assert(haplotype);
    assert(refHaplotype);
    assert( haplotype->getCigar());
    assert(aligner);
    assert(referenceStart >= 0);

    // compute the smith-waterman alignment of read -> haplotype
    auto swPairwiseAlignment = aligner->align(haplotype->getBases(), haplotype->getLength(), originalRead->getBasesNoCopy(), originalRead->getLength(), const_cast<SWParameters *>(&CigarUtils::ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS), SWOverhangStrategy::SOFTCLIP);

    if ( swPairwiseAlignment->getAlignmentOffset() == -1 ) {
        // sw can fail (reasons not clear) so if it happens just don't realign the read
        return originalRead;
    }

    auto swCigar = consolidateCigar(swPairwiseAlignment->getCigar());

    // since we're modifying the read we need to clone it
    std::shared_ptr<SAMRecord> read = std::make_shared<SAMRecord>(*originalRead);
    // only informative reads are given the haplotype tag to enhance visualization
//    if ( isInformative ) {
//        read.setAttribute(HAPLOTYPE_TAG, haplotype.hashCode());
//    }

    // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
    auto extendedHaplotypeCigar = haplotype->getConsolidatedPaddedCigar(1000);
    int readStartOnHaplotype = calcFirstBaseMatchingReferenceInCigar(extendedHaplotypeCigar, swPairwiseAlignment->getAlignmentOffset());
    int readStartOnReference = referenceStart + haplotype->getAlignmentStartHapwrtRef() + readStartOnHaplotype;

    // compute the read -> ref alignment by mapping read -> hap -> ref from the
    // SW of read -> hap mapped through the given by hap -> ref
    auto haplotypeToRef = trimCigarByBases(extendedHaplotypeCigar, swPairwiseAlignment->getAlignmentOffset(), extendedHaplotypeCigar->getReadLength() - 1);
    auto readToRefCigarRaw = applyCigarToCigar(swCigar, haplotypeToRef);
    auto readToRefCigarClean = cleanUpCigar(readToRefCigarRaw);
    auto readToRefCigar =  leftAlignIndel(readToRefCigarClean, refHaplotype->getBases(), refHaplotype->getLength(),
                                          originalRead->getBasesNoCopy(), originalRead->getLength(), readStartOnHaplotype, 0, true);

    int leadingDeletions = readToRefCigarClean->getReferenceLength() - readToRefCigar->getReferenceLength();
    read->setPosition(read->getContig(), readStartOnReference + leadingDeletions);

    // the SW Cigar does not contain the hard clips of the original read
    auto& originalCigar = originalRead->getCigar();
    CigarElement firstElement = originalCigar->getFirstCigarElement();
    CigarElement lastElement = originalCigar->getLastCigarElement();
    std::vector<CigarElement> readToRefCigarElementsWithHardClips;
    if(firstElement.getOperator() == CigarOperator::H)
        readToRefCigarElementsWithHardClips.push_back(firstElement);
    for(auto & cigarElement: readToRefCigar->getCigarElements())
        readToRefCigarElementsWithHardClips.emplace_back(cigarElement);
    if(lastElement.getOperator() == CigarOperator::H)
        readToRefCigarElementsWithHardClips.push_back(lastElement);

    read->setCigar(std::make_shared<Cigar>(readToRefCigarElementsWithHardClips));
    if(readToRefCigar->getReadLength() != read->getLength())
        throw std::exception();

    delete swPairwiseAlignment;
    return read;
}

int
AlignmentUtils::calcFirstBaseMatchingReferenceInCigar(const std::shared_ptr<Cigar> &cigar, int readStartByBaseOfCigar)
{
    assert(cigar != nullptr);
    assert(readStartByBaseOfCigar <= cigar->getReadLength());

    int hapOffset = 0, refOffset = 0;
    for (CigarElement& ce : cigar->getCigarElements() ) {
        for ( int i = 0; i < ce.getLength(); i++ ) {
            switch ( ce.getOperator() ) {
                case M:
                case EQ:
                case X:
                    if ( hapOffset >= readStartByBaseOfCigar )
                        return refOffset;
                    hapOffset++;
                    refOffset++;
                    break;
                case I:
                case S:
                    hapOffset++;
                    break;
                case D:
                    refOffset++;
                    break;
                default:
                    throw "calcFirstBaseMatchingReferenceInCigar does not support cigar";
            }
        }
    }
    throw std::exception();
}

std::shared_ptr<Cigar> AlignmentUtils::applyCigarToCigar(const std::shared_ptr<Cigar> &firstToSecond,
                                                                const std::shared_ptr<Cigar> &secondToThird) {
    std::vector<CigarElement> newElements;
    int nElements12 = firstToSecond->numCigarElements();
    int nElements23 = secondToThird->numCigarElements();

    int cigar12I = 0, cigar23I = 0;
    int elt12I = 0, elt23I = 0;

    while ( cigar12I < nElements12 && cigar23I < nElements23 ) {
        auto & elt12 = firstToSecond->getCigarElement(cigar12I);
        auto & elt23 = secondToThird->getCigarElement(cigar23I);

        CigarPairTransform transform = getTransformer(elt12.getOperator(), elt23.getOperator());

        if ( transform.op13 != CigarOperator::_NULL ) // skip no ops
            newElements.emplace_back(CigarElement(1, transform.op13));

        elt12I += transform.advance12;
        elt23I += transform.advance23;

        // if have exhausted our current element, advance to the next one
        if ( elt12I == elt12.getLength() ) { cigar12I++; elt12I = 0; }
        if ( elt23I == elt23.getLength() ) { cigar23I++; elt23I = 0; }
    }
    return AlignmentUtils::consolidateCigar(std::make_shared<Cigar>(newElements));
}

CigarPairTransform AlignmentUtils::getTransformer(CigarOperator op12, CigarOperator op23)
{
    for(CigarPairTransform& transform : cigarPairTransformers)
    {
        if ( transform.op12.find(op12) != transform.op12.end() && transform.op23.find(op23) != transform.op23.end())
            return transform;
    }
    throw std::exception();
}

bool AlignmentUtils::isInsideDeletion(std::shared_ptr<Cigar> cigar, int offset) {
    assert(cigar != nullptr);
    if(offset < 0)
        return false;

    // pos counts read bases
    int pos = 0;
    int prevPos = 0;

    for (CigarElement& ce : cigar->getCigarElements()) {

        switch (ce.getOperator()) {
            case I:
            case S:
            case D:
            case M:
            case EQ:
            case X:
                prevPos = pos;
                pos += ce.getLength();
                break;
            case H:
            case P:
            case N:
                break;
            default:
                throw "Unsupported cigar operator: ";
        }

        // Is the offset inside a deletion?
        if ( prevPos < offset && pos >= offset && ce.getOperator() == CigarOperator::D ) {
            return true;
        }
    }
    return false;
}

int AlignmentUtils::calcAlignmentByteArrayOffset(std::shared_ptr<Cigar> cigar, int offset, bool isDeletion,
                                                 int alignmentStart, int refLocus) {
    if ( cigar == nullptr ) throw std::invalid_argument("attempting to find the alignment position from a CIGAR that is null");
    if ( offset < -1 ) throw std::invalid_argument("attempting to find the alignment position with an offset that is negative (and not -1)");
    if ( alignmentStart < 0 ) throw std::invalid_argument("attempting to find the alignment position from an alignment start that is negative");
    if ( refLocus < 0 ) throw std::invalid_argument("attempting to find the alignment position from a reference position that is negative");
    if ( offset >= cigar->getReadLength() ) throw std::invalid_argument("attempting to find the alignment position of an offset than is larger than the read length");

    int pileupOffset = offset;

    // Reassign the offset if we are in the middle of a deletion because of the modified representation of the read bases
    if (isDeletion) {
        pileupOffset = refLocus - alignmentStart;
        CigarElement& ce = cigar->getCigarElement(0);
        if (ce.getOperator() == CigarOperator::S) {
            pileupOffset += ce.getLength();
        }
    }

    int pos = 0;
    int alignmentPos = 0;

    for (int iii = 0; iii < cigar->numCigarElements(); iii++) {
        CigarElement& ce = cigar->getCigarElement(iii);
        int elementLength = ce.getLength();

        switch (ce.getOperator()) {
            case I:
            case S: // TODO -- I don't think that soft clips should be treated the same as inserted bases here. Investigation needed.
                pos += elementLength;
                if (pos >= pileupOffset) {
                    return alignmentPos;
                }
                break;
            case D:
                if (!isDeletion) {
                    alignmentPos += elementLength;
                } else {
                    if (pos + elementLength - 1 >= pileupOffset) {
                        return alignmentPos + (pileupOffset - pos);
                    } else {
                        pos += elementLength;
                        alignmentPos += elementLength;
                    }
                }
                break;
            case M:
            case EQ:
            case X:
                if (pos + elementLength - 1 >= pileupOffset) {
                    return alignmentPos + (pileupOffset - pos);
                } else {
                    pos += elementLength;
                    alignmentPos += elementLength;
                }
                break;
                case H:
                    case P:
                        case N:
                            break;
                default:
                    throw "Unsupported cigar operator: ";
        }

    }
    return alignmentPos;
}

int AlignmentUtils::getNumAlignedBasesCountingSoftClips(std::shared_ptr<SAMRecord> r) {
    int n = 0;
    auto& cigar = r->getCigar();
    if (cigar == nullptr)
        return 0;

    for (CigarElement& e : cigar->getCigarElements())
        if (ALIGNED_TO_GENOME_PLUS_SOFTCLIPS.find(e.getOperator()) != ALIGNED_TO_GENOME_PLUS_SOFTCLIPS.end())
            n += e.getLength();

    return n;
}