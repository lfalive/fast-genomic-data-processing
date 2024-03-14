//
// Created by 梦想家xixi on 2021/11/25.
//

#include <cstring>
#include <algorithm>
#include <memory>
#include <htslib/sam.h>
#include "CigarUtils.h"
#include "Mutect2Utils.h"
#include "SWNativeAlignerWrapper.h"
#include "AlignmentUtils.h"

const SWParameters  CigarUtils::NEW_SW_PARAMETERS = SWParameters(200, -150, -260, -11);
const SWParameters  CigarUtils::ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS = SWParameters(10, -15, -30, -5);

std::shared_ptr<Cigar> CigarUtils::calculateCigar(const std::shared_ptr<uint8_t[]>& refSeq, int refLength, const std::shared_ptr<uint8_t[]>& altSeq, int altLength) {
    Mutect2Utils::validateArg(refSeq.get(), "refSeq");
    Mutect2Utils::validateArg(altSeq.get(), "altSeq");
    if(altLength == 0) {
        std::vector<CigarElement> elements;
        elements.emplace_back(CigarElement(refLength, D));
        return std::make_shared<Cigar>(elements);
    }
    if(altLength == refLength) {
        int mismatchCount = 0;
        uint8_t * alt = altSeq.get();
        uint8_t * ref = refSeq.get();
        for (int n = 0; n < refLength && mismatchCount <= 2; n++) {
            mismatchCount += (alt[n] == ref[n] ? 0 : 1);
        }
        if(mismatchCount <= 2) {
            std::shared_ptr<Cigar> matching(new Cigar());
            matching->add(CigarElement(refLength, M));
            return matching;
        }
    }
    std::shared_ptr<Cigar> nonStandard;
    int paddedRefLength = refLength + 2*SW_PAD;
    std::shared_ptr<uint8_t[]> paddedRef(new uint8_t[paddedRefLength]);
    int paddedPathLength = altLength + 2*SW_PAD;
    std::shared_ptr<uint8_t[]> paddedPath(new uint8_t[paddedPathLength]);
    uint8_t tmp[10] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};
    memcpy(paddedRef.get(), tmp, 10);
    memcpy(paddedRef.get()+10, refSeq.get(), refLength);
    memcpy(paddedRef.get()+10+refLength, tmp, 10);
    memcpy(paddedPath.get(), tmp, 10);
    memcpy(paddedPath.get()+10, altSeq.get(), altLength);
    memcpy(paddedPath.get()+10+altLength, tmp, 10);
    SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
    SmithWatermanAlignment* alignment = wrapper.align(paddedRef, paddedRefLength, paddedPath, paddedPathLength,
                                                      const_cast<SWParameters *>(&NEW_SW_PARAMETERS), SOFTCLIP);
    if ( isSWFailure(alignment) ) {
		delete alignment;
        return nullptr;
    }

    int baseStart = SW_PAD;
    int baseEnd = paddedPathLength - SW_PAD - 1;
    nonStandard = AlignmentUtils::trimCigarByBases(alignment->getCigar(), baseStart, baseEnd);
    if(nonStandard->getReferenceLength() != refLength) {
        nonStandard->add(CigarElement(refLength - nonStandard->getReferenceLength(), D));
    }

	delete alignment;
    return leftAlignCigarSequentially(nonStandard, refSeq, refLength, altSeq, altLength, 0, 0);
}

bool CigarUtils::isSWFailure(SmithWatermanAlignment *alignment) {
    if(alignment->getAlignmentOffset() > 0) {
        return true;
    }
    for(const CigarElement & ce : alignment->getCigar()->getCigarElements()) {
        if(ce.getOperator() == S)
            return true;
    }
    //std::any_of(alignment->getCigar()->getCigarElements().begin(), alignment->getCigar()->getCigarElements().end(), [] (CigarElement& ce) {ce.getOperator() == S; return true;});
    return false;
}

std::shared_ptr<Cigar>
CigarUtils::leftAlignCigarSequentially(std::shared_ptr<Cigar> & cigar, const std::shared_ptr<uint8_t[]>& refSeq, int refLength, const std::shared_ptr<uint8_t[]>& readSeq, int readLength,
                                       int refIndex, int readIndex) {
    Mutect2Utils::validateArg(cigar.get(), "cigar null");
    Mutect2Utils::validateArg(refSeq.get(), "refSeq null");
    Mutect2Utils::validateArg(readSeq.get(), "readSeq null");

    std::shared_ptr<Cigar> cigarToReturn(new Cigar());
    std::shared_ptr<Cigar> cigarToAlign(new Cigar());

    for(int i = 0; i < cigar->numCigarElements(); i++) {
        CigarElement ce = cigar->getCigarElement(i);
        if(ce.getOperator() == D || ce.getOperator() == I) {
            cigarToAlign->add(ce);
            std::shared_ptr<Cigar> leftAligned = AlignmentUtils::leftAlignSingleIndel(cigarToAlign, refSeq, refLength, readSeq, readLength, refIndex, readIndex,
                                                                      false);
            for(CigarElement toAdd : leftAligned->getCigarElements()) {cigarToReturn->add(toAdd);}
            refIndex += cigarToAlign->getReferenceLength();
            readIndex += cigarToAlign->getReadLength();
            cigarToAlign = std::make_shared<Cigar>();
        } else {
            cigarToAlign->add(ce);
        }
    }
    if(!cigarToAlign->isEmpty()) {
        for(CigarElement toAdd : cigarToAlign->getCigarElements()) {
            cigarToReturn->add(toAdd);
        }
    }

    std::shared_ptr<Cigar> result = AlignmentUtils::consolidateCigar(cigarToReturn);
    Mutect2Utils::validateArg(result->getReferenceLength() == cigar->getReferenceLength(), "leftAlignCigarSequentially failed to produce a valid CIGAR.");
    return result;
}

bool CigarUtils::isGood(const std::shared_ptr<Cigar> & c) {
    Mutect2Utils::validateArg(c.get(), "cigar is null");

    if(!c->isValid("", -1).empty()) {
        return false;
    }
    std::vector<CigarElement>& elems = c->getCigarElements();
    if(hasConsecutiveIndels(elems)) {
        return false;
    }
    if(startsWithDeletionIgnoringClips(elems)) {
        return false;
    }
    std::vector<CigarElement> elemsRev;
    unsigned size = elems.size();
    elemsRev.reserve(size);
    for(unsigned i = 0; i < size; i++) {
        elemsRev.emplace_back(elems[size -1 -i]);
    }
    return !startsWithDeletionIgnoringClips(elemsRev);
}

bool CigarUtils::isGood(bam1_t * read)
{
    if(!Cigar::isValid(read)){
        return false;
    }
    uint32_t * elems = bam_get_cigar(read);
    uint32_t n_cigar = read->core.n_cigar;
    if(hasConsecutiveIndels(elems, n_cigar)){
        return false;
    }
    if(startsWithDeletionIgnoringClips(elems, n_cigar)){
        return false;
    }

    uint32_t elemsRev[n_cigar];
    for(uint32_t i=0; i<n_cigar; i++){
        elemsRev[i] = elems[n_cigar - i - 1];
    }
    return !startsWithDeletionIgnoringClips(elemsRev, n_cigar);
}

bool CigarUtils::hasConsecutiveIndels(std::vector<CigarElement> &elems) {
    bool prevIndel = false;
    for(CigarElement elem : elems) {
        CigarOperator op = elem.getOperator();
        bool isIndel = (op == I || op == D);
        if(prevIndel && isIndel) {
            return true;
        }
        prevIndel = isIndel;
    }
    return false;
}

bool CigarUtils::hasConsecutiveIndels(uint32_t * cigarArray, uint32_t n_cigar){
    bool prevIndel = false;
    for(uint32_t i=0; i<n_cigar; i++)
    {
        uint32_t op = bam_cigar_op(cigarArray[i]);
        bool isIndel = (op == BAM_CINS || op == BAM_CDEL);
        if(prevIndel && isIndel) {
            return true;
        }
        prevIndel = isIndel;
    }
    return false;
}

bool CigarUtils::startsWithDeletionIgnoringClips(uint32_t * cigarArray, uint32_t n_cigar){
    bool isClip = true;
    int i = 0;
    uint32_t op;
    while(i < n_cigar && isClip) {
        uint32_t elem = cigarArray[i];
        op = bam_cigar_op(elem);
        isClip = (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP);
        i++;
    }
    return op == BAM_CDEL;
}

bool CigarUtils::startsWithDeletionIgnoringClips(const std::vector<CigarElement> &elems) {
    bool isClip = true;
    CigarOperator op;
    int i = 0;
    while(i < elems.size() && isClip) {
        CigarElement elem = elems[i];
        op = elem.getOperator();
        isClip = (op == H || op == S);
        i++;
    }
    return op == D;
}

bool CigarUtils::containsNOperator(std::vector<CigarElement> cigarElements) {
    return std::any_of(cigarElements.begin(), cigarElements.end(), [](CigarElement & cigarElement)->bool { return cigarElement.getOperator() == N;});
}

bool CigarUtils::containsNOperator(uint32_t n_cigar, uint32_t* cigarArray){
    for(uint32_t i=0; i<n_cigar; i++)
    {
        if(bam_cigar_op(cigarArray[i]) == BAM_CREF_SKIP)
            return true;
    }
    return false;
}
