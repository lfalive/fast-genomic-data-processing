//
// Created by 梦想家xixi on 2021/12/18.
//

#include <cmath>
#include <cassert>
#include "ReadUtils.h"
#include "samtools/SAMUtils.h"
#include "ReadConstants.h"

std::string ReadUtils::BQSR_BASE_DELETION_QUALITIES = "BD";
std::string ReadUtils::BQSR_BASE_INSERTION_QUALITIES = "BI";
const int ReadUtils::CLIPPING_GOAL_NOT_REACHED = -1;

std::shared_ptr<SAMRecord> ReadUtils::emptyRead(SAMRecord &read)
{
    std::shared_ptr<SAMRecord> emptyRead(new SAMRecord(read));
    emptyRead->setIsUnmapped();
    emptyRead->setMappingQuality(0);
    std::shared_ptr<Cigar> newCigar(new Cigar());
    emptyRead->setCigar(newCigar);
    emptyRead->setBases(nullptr, 0);
    emptyRead->setBaseQualities(nullptr, 0);
    emptyRead->clearAttributes();

    return emptyRead;
}

std::shared_ptr<SAMRecord> ReadUtils::emptyRead(std::shared_ptr<SAMRecord> & read) {
    std::shared_ptr<SAMRecord> emptyRead(new SAMRecord(*read));
    emptyRead->setIsUnmapped();
    emptyRead->setMappingQuality(0);
    std::shared_ptr<Cigar> newCigar(new Cigar());
    emptyRead->setCigar(newCigar);
    emptyRead->setBases(nullptr, 0);
    emptyRead->setBaseQualities(nullptr, 0);
    emptyRead->clearAttributes();

//    std::string readGroup = read->getReadGroup();
//    if(readGroup.empty()) {
//        emptyRead->setAttribute((std::string&)"RG", readGroup);
//    }
    return emptyRead;
}

void ReadUtils::assertAttributeNameIsLegal(std::string &attributeName) {
    if(attributeName.empty() || attributeName.length() !=2) {
        throw std::invalid_argument("Read attribute invalid: attribute names must be non-null two-character Strings matching the pattern /[A-Za-z][A-Za-z0-9]/");
    }
}

int
ReadUtils::getReadCoordinateForReferenceCoordinate(int alignmentStart, std::shared_ptr<Cigar> cigar, int refCoord, ClippingTail tail,
                                                   bool allowGoalNotReached) {
    std::pair<int, bool> result = getReadCoordinateForReferenceCoordinate(alignmentStart, cigar, refCoord, allowGoalNotReached);
    int readCoord = result.first;
    if(result.second && tail == RIGHT_TAIL){
        readCoord++;
    }
    CigarElement* firstElementIsInsertion = readStartsWithInsertion(cigar);
    if(readCoord == 0 && tail == LEFT_TAIL && firstElementIsInsertion != nullptr) {
        readCoord = std::min(firstElementIsInsertion->getLength(), cigar->getReadLength() - 1);
    }
    return readCoord;
}

std::pair<int, bool> ReadUtils::getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> &read, int refCoord)
{
    return getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), refCoord, false);
}

std::pair<int, bool> ReadUtils::getReadCoordinateForReferenceCoordinate(int alignmentStart, std::shared_ptr<Cigar> cigar, int refCoord,
                                                                        bool allowGoalNotReached) {
    int readBases = 0;
    int refBases = 0;
    bool fallsInsideDeletionOrSkippedRegion = false;
    bool endJustBeforeDeletionOrSkippedRegion = false;
    bool fallsInsideOrJustBeforeDeletionOrSkippedRegion = false;
    int goal = refCoord - alignmentStart;

    if (goal < 0) {
        if (allowGoalNotReached) {
            return {CLIPPING_GOAL_NOT_REACHED, false};
        } else {
            throw std::invalid_argument("Somehow the requested coordinate is not covered by the read. Too many deletions?");
        }
    }
    bool goalReached = refBases == goal;
    std::vector<CigarElement> cigars = cigar->getCigarElements();
    for(int i = 0; !goalReached && i < cigars.size();) {
        CigarElement cigarElement = cigars[i];
        i++;
        int shift = 0;
        if(CigarOperatorUtils::getConsumesReferenceBases(cigarElement.getOperator()) || cigarElement.getOperator() == S) {
            if (refBases + cigarElement.getLength() < goal) {
                shift = cigarElement.getLength();
            } else {
                shift = goal - refBases;
            }

            refBases += shift;
        }
        goalReached = refBases == goal;
        if(!goalReached && CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
            readBases += cigarElement.getLength();
        }
        if(goalReached) {
            bool endsWithinCigar = shift < cigarElement.getLength();
            if(!endsWithinCigar && i >= cigars.size()) {
                if(allowGoalNotReached) {
                    return {CLIPPING_GOAL_NOT_REACHED, false};
                } else {
                    throw std::invalid_argument("Reference coordinate corresponds to a non-existent base in the read.");
                }
            }
            CigarElement* nextCigarElement = nullptr;
            if(endsWithinCigar) {
                fallsInsideDeletionOrSkippedRegion = (cigarElement.getOperator() == D || cigarElement.getOperator() == N);
            } else {
                nextCigarElement = &cigars[i];
                i++;
                if(nextCigarElement->getOperator() == I) {
                    readBases += nextCigarElement->getLength();
                    if(i >= cigars.size()) {
                        if(allowGoalNotReached) {
                            return {CLIPPING_GOAL_NOT_REACHED, false};
                        } else {
                            throw std::invalid_argument("Reference coordinate corresponds to a non-existent base in the read.");
                        }
                    }
                    nextCigarElement = &cigars[i];
                    i++;
                }
                endJustBeforeDeletionOrSkippedRegion = (nextCigarElement->getOperator() == D || nextCigarElement->getOperator() == N);
            }
            fallsInsideOrJustBeforeDeletionOrSkippedRegion = endJustBeforeDeletionOrSkippedRegion || fallsInsideDeletionOrSkippedRegion;
            if(!fallsInsideOrJustBeforeDeletionOrSkippedRegion && CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                readBases += shift;
            } else if(endJustBeforeDeletionOrSkippedRegion && CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                readBases += shift - 1;
            } else if(fallsInsideDeletionOrSkippedRegion ||
                    (endJustBeforeDeletionOrSkippedRegion && nextCigarElement->getOperator() == N) ||
                    (endJustBeforeDeletionOrSkippedRegion && nextCigarElement->getOperator() == D)) {
                readBases--;
            }
        }
    }
    if(!goalReached) {
        if(allowGoalNotReached) {
            return {CLIPPING_GOAL_NOT_REACHED, false};
        } else {
            throw std::invalid_argument("Somehow the requested coordinate is not covered by the read. Alignment");
        }
    }

    return {readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion};
}

CigarElement* ReadUtils::readStartsWithInsertion(const std::shared_ptr<Cigar>& cigarForRead, bool ignoreSoftClipOps) {
    for(CigarElement& cigarElement : cigarForRead->getCigarElements()) {
        if(cigarElement.getOperator() == I) {
            return &cigarElement;
        } else if (cigarElement.getOperator() != H && (!ignoreSoftClipOps || cigarElement.getOperator() != S)) {
            break;
        }
    }
    return nullptr;
}

CigarElement* ReadUtils::readStartsWithInsertion(const std::shared_ptr<Cigar>& cigarForRead) {
    return readStartsWithInsertion(cigarForRead, true);
}

int ReadUtils::getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> & read, int refCoord, ClippingTail tail) {
    //int leftmostSafeVariantPosition = std::max(read->getSoftStart(), refCoord);
    return getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), refCoord, tail,
                                                   false);
}

int ReadUtils::getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> & read, int refCoord, ClippingTail tail, bool allowGoalNotReached)
{
    int leftmostSafeVariantPosition = std::max(read->getSoftStart(), refCoord);
    return getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), leftmostSafeVariantPosition, tail, allowGoalNotReached);
}

int ReadUtils::getSoftStart(std::shared_ptr<SAMRecord> & read) {
    Mutect2Utils::validateArg(read != nullptr, "read");

    int softStart = read->getStart();
    for(CigarElement cig : read->getCigarElements()) {
        CigarOperator op = cig.getOperator();
        if(op == S) {
            softStart -= cig.getLength();
        } else if (op != H) {
            break;
        }
    }
    return softStart;
}

int ReadUtils::getSoftEnd(std::shared_ptr<SAMRecord> & read) {
    Mutect2Utils::validateArg(read != nullptr, "read");

    bool foundAlignedBase = false;
    int softEnd = read->getEnd();
    std::vector<CigarElement> cigs = read->getCigarElements();
    for (int i = cigs.size() - 1; i >= 0; --i) {
        CigarElement cig = cigs[i];
        CigarOperator op = cig.getOperator();

        if (op == S){
            softEnd += cig.getLength();
        } else if (op != H) {
            foundAlignedBase = true;
            break;
        }
    }
    if( !foundAlignedBase ) {
        softEnd = read->getEnd();
    }
    return softEnd;
}

bool ReadUtils::hasBaseIndelQualities(std::shared_ptr<SAMRecord> & read) {
    return read->getAttribute(SAMUtils::makeBinaryTag(BQSR_BASE_INSERTION_QUALITIES)) || read->getAttribute(SAMUtils::makeBinaryTag(BQSR_BASE_DELETION_QUALITIES));
}

std::shared_ptr<uint8_t[]> ReadUtils::getExistingBaseInsertionQualities(SAMRecord &read, int &length)
{
    std::string str = read.getAttributeAsString(BQSR_BASE_INSERTION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

std::shared_ptr<uint8_t[]> ReadUtils::getExistingBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    std::string str = read->getAttributeAsString(BQSR_BASE_INSERTION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

std::shared_ptr<uint8_t[]> ReadUtils::getExistingBaseDeletionQualities(SAMRecord &read, int &length)
{
    std::string str = read.getAttributeAsString(BQSR_BASE_DELETION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

std::shared_ptr<uint8_t[]> ReadUtils::getExistingBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    std::string str = read->getAttributeAsString(BQSR_BASE_DELETION_QUALITIES);
    length = str.length();
    return SAMUtils::fastqToPhred(str);
}

std::shared_ptr<uint8_t[]> ReadUtils::getBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    std::shared_ptr<uint8_t[]> quals = getExistingBaseInsertionQualities(read, length);
    if(quals == nullptr) {
        length = read->getBaseQualitiesLength();
        quals = std::shared_ptr<uint8_t[]>(new uint8_t[length]);
        memset(quals.get(), DEFAULT_INSERTION_DELETION_QUAL, length);
    }
    return quals;
}

std::shared_ptr<uint8_t[]> ReadUtils::getBaseInsertionQualities(SAMRecord &read, int &length) {
    std::shared_ptr<uint8_t[]> quals = getExistingBaseInsertionQualities(read, length);
    if(quals == nullptr) {
        length = read.getBaseQualitiesLength();
        quals = std::shared_ptr<uint8_t[]>(new uint8_t[length]);
        memset(quals.get(), DEFAULT_INSERTION_DELETION_QUAL, length);
    }
    return quals;
}

std::shared_ptr<uint8_t[]> ReadUtils::getBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int &length) {
    std::shared_ptr<uint8_t[]> quals = getExistingBaseDeletionQualities(read, length);
    if(quals == nullptr) {
        length = read->getBaseQualitiesLength();
        quals = std::shared_ptr<uint8_t[]>(new uint8_t[length]);
        memset(quals.get(), DEFAULT_INSERTION_DELETION_QUAL, length);
    }
    return quals;
}

std::shared_ptr<uint8_t[]> ReadUtils::getBaseDeletionQualities(SAMRecord &read, int &length){
    std::shared_ptr<uint8_t[]> quals = getExistingBaseDeletionQualities(read, length);
    if(quals == nullptr) {
        length = read.getBaseQualitiesLength();
        quals = std::shared_ptr<uint8_t[]>(new uint8_t[length]);
        memset(quals.get(), DEFAULT_INSERTION_DELETION_QUAL, length);
    }
    return quals;
}

void ReadUtils::setInsertionBaseQualities(std::shared_ptr<SAMRecord> & read, std::shared_ptr<uint8_t[]> quals, int length) {
    read->setAttribute(BQSR_BASE_INSERTION_QUALITIES, quals == nullptr ? "" : SAMUtils::phredToFastq(quals, length));
}

void ReadUtils::setDeletionBaseQualities(std::shared_ptr<SAMRecord> & read, std::shared_ptr<uint8_t[]> quals, int length) {
    read->setAttribute(BQSR_BASE_DELETION_QUALITIES, quals  == nullptr ? "" : SAMUtils::phredToFastq(quals, length));
}

int ReadUtils::getAssignedReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader *header) {
    return header->getSequenceIndex(read->getAssignedContig());
}

int ReadUtils::getSAMFlagsForRead(std::shared_ptr<SAMRecord> & read) {
    int samFlags = 0;

    if(read->isPaired()) {
        samFlags |= 1;
    }
    if(read->isProperlyPaired()) {
        samFlags |= 2;
    }
    if(read->isUnmapped()) {
        samFlags |= 4;
    }
    if(read->isPaired() && read->mateIsUnmapped()) {
        samFlags |= 8;
    }
    if(!read->isUnmapped() && read->isReverseStrand()) {
        samFlags |= 16;
    }
    if(read->isPaired() && ! read->mateIsUnmapped() && read->mateIsReverseStrand()) {
        samFlags |= 32;
    }
    if(read->isFirstOfPair()) {
        samFlags |= 64;
    }
    if(read->isSecondOfPair()) {
        samFlags |= 128;
    }
    if(read->isSecondaryAlignment()) {
        samFlags |= 256;
    }
    if(read->failsVendorQualityCheck()) {
        samFlags |= 512;
    }
    if(read->isDuplicate()) {
        samFlags |= 1024;
    }
    if(read->isSupplementaryAlignment()) {
        samFlags |= 2048;
    }
    return samFlags;
}

int ReadUtils::getMateReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader *header) {
    if(read->mateIsUnmapped()) {
        return -1;
    }
    return header->getSequenceIndex(read->getMateContig());
}

bool ReadUtils::alignmentAgreesWithHeader(SAMFileHeader *header, std::shared_ptr<SAMRecord> & read) {
    int referenceIndex = getReferenceIndex(read, header);

    if(! read->isUnmapped() && referenceIndex == SAMRecord::NO_ALIGNMENT_REFERENCE_INDEX) {
        return false;
    }
    SAMSequenceRecord contigHeader = header->getSequenceDictionary().getSequences()[referenceIndex];
    return read->isUnmapped() || read->getStart() <= contigHeader.getSequenceLength();
}

bool ReadUtils::alignmentAgreesWithHeader(sam_hdr_t * hdr, bam1_t * read){
    if(!isUnmapped(read, hdr) && read->core.tid == SAMRecord::NO_ALIGNMENT_REFERENCE_INDEX){
        return false;
    }
    return isUnmapped(read, hdr) || read->core.pos <= sam_hdr_tid2len(hdr, read->core.tid);
}

int ReadUtils::getReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader *header) {
//    if(read->isUnmapped()) {
//        return SAMRecord::NO_ALIGNMENT_REFERENCE_INDEX;
//    }
    return header->getSequenceIndex(read->getContig());
}

bool ReadUtils::hasWellDefinedFragmentSize(std::shared_ptr<SAMRecord> & read) {
    if(read->getFragmentLength() == 0) {
        return false;
    }
    if( ! read->isPaired()) {
        return false;
    }
    if(read->isUnmapped() || read->mateIsUnmapped()) {
        return false;
    }
    if(read->isReverseStrand() == read->mateIsReverseStrand()) {
        return false;
    }
    if(read->isReverseStrand()) {
        return read->getEnd() > read->getMateStart();
    } else {
        return read->getStart() <= read->getMateStart() + read->getFragmentLength();
    }
}

bool ReadUtils::hasWellDefinedFragmentSize(SAMRecord * read) {
    if(read->getFragmentLength() == 0) {
        return false;
    }
    if( ! read->isPaired()) {
        return false;
    }
    if(read->isUnmapped() || read->mateIsUnmapped()) {
        return false;
    }
    if(read->isReverseStrand() == read->mateIsReverseStrand()) {
        return false;
    }
    if(read->isReverseStrand()) {
        return read->getEnd() > read->getMateStart();
    } else {
        return read->getStart() <= read->getMateStart() + read->getFragmentLength();
    }
}

bool ReadUtils::hasWellDefinedFragmentSize(bam1_t * read, sam_hdr_t * hdr) {
    if(read->core.isize == 0)
        return false;
    if(!isPaired(read))
        return false;
    if(isUnmapped(read, hdr) || mateIsUnmapped(read, hdr))
        return false;
    if(isReverseStrand(read) == mateIsReverseStrand(read))
        return false;
    if(isReverseStrand(read))
        return getEnd(read) > read->core.mpos;
    else
        return read->core.pos <= read->core.mpos + read->core.isize;
}

int ReadUtils::getAdaptorBoundary(SAMRecord * read) {
    if(!hasWellDefinedFragmentSize(read)) {
        return INT32_MIN;
    } else if (read->isReverseStrand()) {
        return read->getMateStart() - 1;
    } else {
        int insertSize = std::abs(read->getFragmentLength());
        return read->getStart() + insertSize;
    }
}

int ReadUtils::getAdaptorBoundary(bam1_t * read, sam_hdr_t * hdr) {
    if(!hasWellDefinedFragmentSize(read, hdr)){
        return INT32_MIN;
    } else if(isReverseStrand(read)){
        return read->core.mpos - 1;
    } else {
        int insertSize = std::abs(read->core.isize);
        return read->core.pos + insertSize;
    }
}

hts_pos_t ReadUtils::getEnd(bam1_t *read)
{
    return bam_endpos(read) - 1;
}

bool ReadUtils::isPaired(bam1_t *read)
{
    return read->core.flag & 1;
}

bool ReadUtils::getProperPairFlagUnchecked(bam1_t *read)
{
    return read->core.flag & 2;
}

bool ReadUtils::isReverseStrand(bam1_t *read)
{
    return read->core.flag & 16;
}

bool ReadUtils::getMateUnmappedFlagUnchecked(bam1_t *read)
{
    return (read->core.flag & 8) != 0;
}

bool ReadUtils::getMateNegativeStrandFlagUnchecked(bam1_t *read)
{
    return read->core.flag & 32;
}

bool ReadUtils::mateIsUnmapped(bam1_t * read, sam_hdr_t * hdr)
{
    if(!isPaired(read)){
        throw std::invalid_argument("Cannot get mate information for an unpaired read");
    }
    return getMateUnmappedFlagUnchecked(read) || (sam_hdr_tid2name(hdr, read->core.mtid) == nullptr) || strcmp(sam_hdr_tid2name(hdr, read->core.mtid),  SAMRecord::NO_ALIGNMENT_REFERENCE_NAME.c_str()) == 0
        || read->core.mpos == SAMRecord::NO_ALIGNMENT_START;
}

bool ReadUtils::mateIsReverseStrand(bam1_t *read)
{
    return getMateNegativeStrandFlagUnchecked(read);
}

const char * ReadUtils::getReferenceName(bam1_t * read, sam_hdr_t * hdr)
{
    return sam_hdr_tid2name(hdr, read->core.tid);
}

bool ReadUtils::isUnmapped(bam1_t *read, sam_hdr_t * hdr)
{
    const char * refName = getReferenceName(read, hdr);
    return (read->core.flag & BAM_FUNMAP) != 0 ||
    strcmp(refName, SAMRecord::NO_ALIGNMENT_REFERENCE_NAME.c_str()) == 0 ||
    read->core.pos + 1 == SAMRecord::NO_ALIGNMENT_START;
}

bool ReadUtils::isProperlyPaired(bam1_t * read)
{
    return isPaired(read) && getProperPairFlagUnchecked(read);
}

uint8_t ReadUtils::decodeBase(uint8_t base)
{
    switch (base) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            return '-';
    }
}



bool ReadUtils::consumesReadBases(uint32_t cigarElement)
{
    return bam_cigar_type(bam_cigar_op(cigarElement)) & 1;
}

bool ReadUtils::consumesReferenceBases(uint32_t cigarElement)
{
    return bam_cigar_type(bam_cigar_op(cigarElement)) & 2;
}



int ReadUtils::calculateAlignmentStartShift(int n_cigar, uint32_t *oldCigar, int newReadBasesClipped)
{
    int readBasesClipped = 0; // The number of read bases consumed on the new cigar before reference bases are consumed
    int refBasesClipped = 0; // A measure of the reference offset between the oldCigar and the clippedCigar

    bool truncated=false;
    int i=0;

    for (; i<n_cigar; i++)
    {
        uint32_t e = oldCigar[i];

        int curRefLength = bam_cigar_oplen(e);
        int curReadLength = consumesReadBases(e) ? bam_cigar_oplen(e) : 0;

        truncated = readBasesClipped + curReadLength > newReadBasesClipped;
        if (truncated) {
            curReadLength = newReadBasesClipped - readBasesClipped;
            curRefLength = curReadLength;
        }

        if (!consumesReferenceBases(e))
            curRefLength = 0;

        readBasesClipped += curReadLength;
        refBasesClipped += curRefLength;

        if (readBasesClipped >= newReadBasesClipped || truncated) {
            break;
        }
    }

    // needed only if the clipping ended at a cigar element boundary and is followed by either N or D
    if (readBasesClipped == newReadBasesClipped && !truncated) {
        while (i < n_cigar - 1){
            i++;

            if (consumesReadBases(oldCigar[i]) || !consumesReferenceBases(oldCigar[i]))
                break;

            refBasesClipped += bam_cigar_oplen(oldCigar[i]);
        }
    }

    return refBasesClipped;
}

bool ReadUtils::isBaseInsideAdaptor(std::shared_ptr<SAMRecord> & read, long basePos) {
    int adaptorBoundary = read->getAdaptorBoundary();
    if(adaptorBoundary == INT32_MIN || read->getFragmentLength() > 100)
        return false;
    return read->isReverseStrand() ? basePos <= adaptorBoundary : basePos >= adaptorBoundary;
}

bool ReadUtils::isInsideRead(std::shared_ptr<SAMRecord> &read, int referenceCoordinate) {
    return referenceCoordinate >= read->getStart() && referenceCoordinate <= read->getEnd();
}

bool ReadUtils::isF2R1(std::shared_ptr<SAMRecord> read) {
    return read->isReverseStrand() == read->isFirstOfPair();
}

