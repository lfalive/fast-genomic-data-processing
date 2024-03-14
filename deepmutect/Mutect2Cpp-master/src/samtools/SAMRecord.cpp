//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMRecord.h"
#include "SAMFlag.h"
#include "ReadConstants.h"
#include "SAMUtils.h"
#include "ReadUtils.h"
#include <memory>
#include <iostream>
#include <utility>

const std::string SAMRecord::NO_ALIGNMENT_REFERENCE_NAME = "*";

std::string &SAMRecord::getName() {
    return mReadName;
}

void SAMRecord::setName(std::string& name) {
    mReadName = name;
}

int SAMRecord::getLength() const {
    return baseLength;
}

void SAMRecord::setPosition(std::string& contig, int start) {
    if(contig.empty() || contig == NO_ALIGNMENT_REFERENCE_NAME || start < 0) {
        throw std::invalid_argument("contig must be non-null and start must be >= 0");
    }
    mReferenceName = contig;
    mAlignmentStart = start;
    mAlignmentEnd = start + mCigar->getReferenceLength() - 1;
    setFlag(false, 4);
}

void SAMRecord::setFlag(bool flag, int bit) {
    if(flag) {
        mFlags |= bit;
    } else {
        mFlags &= ~bit;
    }
}

void SAMRecord::setPosition(Locatable *locatable) {
    Mutect2Utils::validateArg(locatable, "Cannot set read position to null");
    std::string contig = locatable->getContig();
    setPosition(contig, locatable->getStart());
}

std::string &SAMRecord::getAssignedContig() {
    return mReferenceName;
}

int SAMRecord::getAssignedStart() const {
    return mAlignmentStart;
}

bool SAMRecord::getReadUnmappedFlag() const {
    return (mFlags & 4) != 0;
}

bool SAMRecord::isUnmapped() {
    return getReadUnmappedFlag() || mReferenceName.empty() || mReferenceName == NO_ALIGNMENT_REFERENCE_NAME ||
    mAlignmentStart + 1 == NO_ALIGNMENT_START;
}

int SAMRecord::getUnclippedStart() {
    if(isUnmapped()) {
        return ReadConstants::UNSET_POSITION;
    }
    return SAMUtils::getUnclippedStart(mAlignmentStart, mCigar);
}

int SAMRecord::getUnclippedEnd() {
    if(isUnmapped())
        return ReadConstants::UNSET_POSITION;

    return SAMUtils::getUnclippedEnd(mAlignmentEnd, mCigar);
}

std::string& SAMRecord::getMateContig() {
    if(isUnmapped())
        return (std::string&)"";
    return mMateReferenceName;
}

bool SAMRecord::isPaired() const {
    return (mFlags & 1) != 0;
}

int SAMRecord::getMateStart() {
    if(mateIsUnmapped()) {
        return ReadConstants::UNSET_POSITION;
    }
    return mMateAlignmentStart;
}

bool SAMRecord::mateIsUnmapped() {
    if(!isPaired()) {
        throw std::invalid_argument("Cannot get mate information for an unpaired read");
    }

    return getMateUnmappedFlag() || mMateReferenceName.empty() || mMateReferenceName == NO_ALIGNMENT_REFERENCE_NAME
    || mMateAlignmentStart + 1 == NO_ALIGNMENT_START;
}

void SAMRecord::requireReadPaired() const {
    if((mFlags & 1) == 0) {
        throw std::invalid_argument("Inappropriate call if not paired read");
    }
}

bool SAMRecord::getMateUnmappedFlagUnchecked() const {
    return (mFlags & 8) != 0;
}

bool SAMRecord::getMateUnmappedFlag() {
    requireReadPaired();
    return getMateUnmappedFlagUnchecked();
}

void SAMRecord::setMatePosition(std::string &contig, int start) {
    if(!contig.empty() || contig == NO_ALIGNMENT_REFERENCE_NAME || start < 0) {
        throw std::invalid_argument("contig must be non-null and start must be >= 0");
    }

    setIsPaired(true);
    //TODO::实现setMateReferenceName
    mMateReferenceName = contig;
    mMateAlignmentStart = start;
    setFlag(false, 8);
}

void SAMRecord::setReadPairedFlag(bool flag) {
    setFlag(flag, 1);
}

void SAMRecord::setProperPairFlag(bool flag) {
    setFlag(flag, 2);
}

void SAMRecord::setIsPaired(bool isPaired) {
    setReadPairedFlag(isPaired);
    if(! isPaired) {
        setProperPairFlag(false);
    }
}

void SAMRecord::setMatePosition(Locatable *locatable) {
    Mutect2Utils::validateArg(locatable, "Cannot set mate position to null");
    std::string contig = locatable->getContig();
    setMatePosition(contig, locatable->getStart());
}

int SAMRecord::getFragmentLength() const {
    return mInferredInsertSize;
}

void SAMRecord::setFragmentLength(int fragmentLength) {
    mInferredInsertSize = fragmentLength;
}

int SAMRecord::getMappingQuality() const {
    return mMappingQuality != NO_MAPPING_QUALITY ? mMappingQuality : NO_MAPPING_QUALITY;
}

void SAMRecord::setMappingQuality(int mappingQuality) {
    if(!(mappingQuality >= 0 && mappingQuality <= 255)) {
        throw std::invalid_argument("mapping quality must be >= 0 and <= 255");
    }
    mMappingQuality = mappingQuality;
}

std::shared_ptr<uint8_t[]>SAMRecord::getBases() {
    if(mReadBases != nullptr){
        std::shared_ptr<uint8_t[]> ret(new uint8_t[baseLength+1]{0});
        std::copy(mReadBases.get(), mReadBases.get() + baseLength, ret.get());
        return ret;
    } else {
        return nullptr;
    }
}

std::shared_ptr<uint8_t[]> & SAMRecord::getBasesNoCopy() {
    return mReadBases;
}

void SAMRecord::setBases(std::shared_ptr<uint8_t[]>bases, int length) {
    mReadBases = std::move(bases);
    baseLength = length;
}

std::shared_ptr<uint8_t[]>SAMRecord::getBaseQualities() {
    if(mBaseQualities != nullptr){
        std::shared_ptr<uint8_t[]> ret{new uint8_t[baseQualitiesLength+1]{0}};
        std::copy(mBaseQualities.get(), mBaseQualities.get()+baseQualitiesLength, ret.get());
        return ret;
    } else {
        return nullptr;
    }
}

std::shared_ptr<uint8_t[]> & SAMRecord::getBaseQualitiesNoCopy() {
    return mBaseQualities;
}

int SAMRecord::getBaseQualitiesLength() const {
    return baseQualitiesLength;
}

void SAMRecord::setBaseQualities(std::shared_ptr<uint8_t[]>baseQualities, int length) {
    mBaseQualities = std::move(baseQualities);
    baseQualitiesLength = length;
}

const std::shared_ptr<Cigar> & SAMRecord::getCigar() {
    //TODO:验证是否需要返回拷贝后的cigar
    return mCigar;
}

std::vector<CigarElement> & SAMRecord::getCigarElements() {
    return mCigar->getCigarElements();
}

CigarElement SAMRecord::getCigarElement(const int index) {
    return mCigar->getCigarElement(index);
}

void SAMRecord::setCigar(std::shared_ptr<Cigar> cigar) {
    //TODO::实现AlignmentBlock
    isCalAdaptorBoundary = false;
    mCigar = std::move(cigar);
    mAlignmentEnd = mAlignmentStart + mCigar->getReferenceLength() - 1;
}

bool SAMRecord::getProperPairFlag() {
    requireReadPaired();
    return getProperPairFlagUnchecked();
}

bool SAMRecord::getProperPairFlagUnchecked() const {
    return (mFlags & 2) != 0;
}

bool SAMRecord::isProperlyPaired() {
    return isPaired() && getProperPairFlag();
}

void SAMRecord::setIsProperlyPaired(bool isProperlyPaired) {
    if(isProperlyPaired) {
        setIsPaired(true);
    }

    setProperPairFlag(isProperlyPaired);
}

SAMRecord::SAMRecord(std::shared_ptr<uint8_t[]>base, int baseLength, std::shared_ptr<uint8_t[]>baseQualities, int baseQualitiesLength,
                     std::string &name) : mReadBases(std::move(base)), baseLength(baseLength), mBaseQualities(std::move(baseQualities)),baseQualitiesLength(baseQualitiesLength), mReadName(name){}

int SAMRecord::getStart() const {
//    if(isUnmapped()) {
//        return ReadConstants::UNSET_POSITION;
//    }
    return mAlignmentStart;
}

int SAMRecord::getEnd() {
//    if(isUnmapped()) {
//        return ReadConstants::UNSET_POSITION;
//    }

    return mAlignmentEnd;
}

void SAMRecord::setIsUnmapped() {
    setReadUnmappedFlag(true);
}

void SAMRecord::setReadUnmappedFlag(bool flag) {
    setFlag(flag, 4);
}

void SAMRecord::clearAttributes() {
    mAttributes = nullptr;
}

void SAMRecord::setAttribute(short tag, void *value, Void_Type type, int length, bool isUnsignedArray) {
    if(value == nullptr) {
        if(mAttributes != nullptr) {
            mAttributes = SAMBinaryTagAndValue::remove(mAttributes, tag);
        }
    } else {
        std::shared_ptr<SAMBinaryTagAndValue> tmp;
        if(!isUnsignedArray) {
            tmp = std::make_shared<SAMBinaryTagAndValue>(tag, value, type, length);
        } else {
            if(type == Short_Array_Type || type == Int_Array_Type || type == Float_Array_Type || type == Uint8_t_Array_Type) {
                throw std::invalid_argument("Attribute type cannot be encoded as an unsigned array.");
            } else {
                tmp = std::make_shared<SAMBinaryTagAndValue>(tag, value, type, length);
            }
        }

        if(mAttributes == nullptr) {
            mAttributes = tmp;
        } else {
            mAttributes = SAMBinaryTagAndValue::insert(mAttributes, tmp);
        }
    }
}

void *SAMRecord::getAttribute(short tag) {
    if(mAttributes == nullptr) {
        return nullptr;
    } else {
        SAMBinaryTagAndValue* tmp = mAttributes->find(tag);
        return tmp != nullptr ? tmp->value : nullptr;
    }
}

std::string &SAMRecord::getReadGroup() {
    return *(std::string*)getAttribute(SAMUtils::makeBinaryTag((std::string &) "RG"));
}

void SAMRecord::setAttribute(std::string &attributeName, const std::string& attributeValue) {
    ReadUtils::assertAttributeNameIsLegal(attributeName);
    //std::cout << attributeName << " " << attributeValue << std::endl;
    setAttribute(attributeName, new std::string(attributeValue), String_Type, 0);
}

void SAMRecord::setAttribute(short tag, void *value, Void_Type type, int length) {
    setAttribute(tag, value, type, length, false);
}

void SAMRecord::setAttribute(std::string &tag, void *value, Void_Type type, int length) {
    setAttribute(SAMUtils::makeBinaryTag(tag), value, type, length);
}

int SAMRecord::getSoftStart() {
    std::shared_ptr<SAMRecord> ptr(new SAMRecord(*this));
    return ReadUtils::getSoftStart(ptr);
}

int SAMRecord::getSoftEnd() {
    std::shared_ptr<SAMRecord> ptr(new SAMRecord(*this));
    return ReadUtils::getSoftStart(ptr);
}

std::string &SAMRecord::getContig() {
    return getReadUnmappedFlag() ? (std::string&)"" : mReferenceName;
}

int SAMRecord::getContigInt() {
	return getReadUnmappedFlag() ? -1 : ContigMap::getContigInt(mReferenceName);
}

std::string SAMRecord::getAttributeAsString(std::string &attributeName) {
    ReadUtils::assertAttributeNameIsLegal(attributeName);
    if(mAttributes == nullptr) {
        return "";
    } else {
        SAMBinaryTagAndValue* tmp = mAttributes->find(SAMUtils::makeBinaryTag(attributeName));
        if(tmp == nullptr)
            return "";
        std::string ret;
        switch (tmp->type) {
            case Uint8_t_Array_Type:
            {
                char* val = (char*) tmp->value;
                if(tmp->length == 0) {
                    ret = "";
                    break;
                }
                char * newVal = new char[tmp->length+1];
                memcpy(newVal, val, tmp->length);
                newVal[tmp->length] = 0;
                ret = std::string(newVal);
                delete[] newVal;
                break;
            }
            case Uint8_Type:
            {
                char* val = (char*) tmp->value;
                ret = *val;
                break;
            }
            case Int_Array_Type:
            {
                std::stringstream ss;
                for(int i = 0; i < tmp->length; i++) {
                    ss << ((int*)tmp->value)[i];
                    ss >> ret;
                }
            }
            case Integer_Type:
            {
                std::stringstream ss;
                ss << *((int*)tmp->value);
                ss >> ret;
            }
            case Float_Array_Type:
            {
                std::stringstream ss;
                for(int i = 0; i < tmp->length; i++) {
                    ss << ((float *)tmp->value)[i];
                    ss >> ret;
                }
            }
            case Float_Type:
            {
                std::stringstream ss;
                ss << *((float *)tmp->value);
                ss >> ret;
            }
            case Long_Type:
            {
                std::stringstream ss;
                ss << *((long *)tmp->value);
                ss >> ret;
            }
            case String_Type:
            {
                return *((std::string*)tmp->value);
            }
            case Short_Type:
            {
                std::stringstream ss;
                ss << *((short *)tmp->value);
                ss >> ret;
            }
            case Short_Array_Type:
            {
                std::stringstream ss;
                for(int i = 0; i < tmp->length; i++) {
                    ss << ((short *)tmp->value)[i];
                    ss >> ret;
                }
            }
        }
        return ret;
    }
}

bool SAMRecord::isReverseStrand() const {
    return (mFlags & 16) != 0;
}

bool SAMRecord::mateIsReverseStrand() {
	if (!isPaired())
		throw std::invalid_argument("Cannot get mate information for an unpaired read");
    requireReadPaired();
    return getMateNegativeStrandFlagUnchecked();
}

bool SAMRecord::getMateNegativeStrandFlagUnchecked() const {
    return (mFlags & 32) != 0;
}

bool SAMRecord::isFirstOfPair() {
    return isPaired() && getFirstOfPairFlag();
}

bool SAMRecord::getFirstOfPairFlag() {
    requireReadPaired();
    return (mFlags & 64) != 0;
}

bool SAMRecord::isSecondOfPair() {
    return isPaired() && getSecondOfPairFlag();
}

bool SAMRecord::getSecondOfPairFlag() {
    requireReadPaired();
    return (mFlags & 128) != 0;
}

bool SAMRecord::isSecondaryAlignment() const {
    return (mFlags & 256) != 0;
}

bool SAMRecord::failsVendorQualityCheck() const {
    return (mFlags & 512) != 0;
}

bool SAMRecord::isDuplicate() const {
    return (mFlags & 1024) != 0;
}

bool SAMRecord::isSupplementaryAlignment() const {
    return (mFlags & 2048) != 0;
}

SAMRecord::SAMRecord(bam1_t *read, sam_hdr_t * hdr, bool load) {
    uint32_t * res = bam_get_cigar(read);
    uint32_t n = read->core.n_cigar;
    std::vector<CigarElement> nCigarElements;
    nCigarElements.reserve(n);
    for(int i = 0; i < n; i++) {
        int length = (int) (res[i] >> 4);
        CigarOperator tmp_cigarOperator = CigarOperatorUtils::binaryToEnum((int) (res[i] & 0xf));
        nCigarElements.emplace_back(CigarElement(length, tmp_cigarOperator));
    }
    mCigar = std::make_shared<Cigar>(nCigarElements);
    mFlags = read->core.flag;
    mMappingQuality = read->core.qual;
    mAlignmentStart = read->core.pos;
    mAlignmentEnd = mAlignmentStart + static_cast<int>(bam_cigar2rlen(n, res)) - 1;
    mReferenceName = std::string(sam_hdr_tid2name(hdr, read->core.tid));
    auto matepointer = sam_hdr_tid2name(hdr, read->core.mtid);
    mMateReferenceName = matepointer != nullptr ? std::string(matepointer) : "";
    mReadName = std::string(bam_get_qname(read));
    baseLength = read->core.l_qseq;
    baseQualitiesLength = read->core.l_qseq;
    if(load) {
        uint8_t * bases = bam_get_seq(read);
        mReadBases = std::shared_ptr<uint8_t[]>(new uint8_t[baseLength+1]{0});
        mBaseQualities = std::shared_ptr<uint8_t[]>(new uint8_t[baseLength+1]{0});
        uint8_t * mReadBases_ = mReadBases.get();
        uint8_t * mBaseQualities_ = mBaseQualities.get();
        for(int i = 0; i < baseLength; i++) {
            uint8_t tmp = bam_seqi(bases, i);
            switch(tmp) {
                case 1 : {
                    mReadBases_[i] = 'A' ;
                    break;
                }
                case 2 : {
                    mReadBases_[i] = 'C';
                    break;
                }
                case 4 : {
                    mReadBases_[i] = 'G';
                    break;
                }
                case 8 : {
                    mReadBases_[i] = 'T';
                    break;
                }
                default : {
                    mReadBases_[i] = 'N';
                    break;
                }
            }
        }
        uint8_t * qual = bam_get_qual(read);
        memcpy(mBaseQualities_, qual, baseQualitiesLength);
    } else {
        mReadBases = nullptr;
        mBaseQualities = nullptr;
    }
    mMateAlignmentStart = static_cast<int>(read->core.mpos);
    mInferredInsertSize = static_cast<int>(read->core.isize);
    mAttributes = nullptr;

    //get the map from pos to cigarElement
    hts_pos_t rlen = bam_cigar2rlen(n, res);
    int start = 0;
    int offset = 0;
    int cigarOffset = -1;
    int currentStart = 0;
    uint32_t cigarLength = 0;
    for(int i=0; i<rlen; i++)
    {
        if(cigarOffset == n-1)
            break;
        while(start <= i)
        {
            cigarOffset++;
            cigarLength = bam_cigar_oplen(res[cigarOffset]);
            if(ReadUtils::consumesReferenceBases(res[cigarOffset]))
            {
                start += cigarLength;
            }
            if(ReadUtils::consumesReadBases(res[cigarOffset]))
            {
                offset += cigarLength;
            }
            currentStart = start - cigarLength;

            PositionToCigarMap.emplace_back(std::pair<int, PositionToCigar>(i, PositionToCigar(cigarOffset, currentStart, offset + i - currentStart - cigarLength)));
        }
    }

}

SAMRecord::~SAMRecord() = default;

int SAMRecord::getAdaptorBoundary() {
    if(isCalAdaptorBoundary)
        return adaptorBoundary;
    adaptorBoundary = ReadUtils::getAdaptorBoundary(this);
    isCalAdaptorBoundary = true;
    return adaptorBoundary;
}

SAMRecord::SAMRecord(const SAMRecord &other) : mFlags(other.mFlags), baseLength(other.baseLength), baseQualitiesLength(other.baseQualitiesLength),
mAlignmentStart(other.mAlignmentStart), mAlignmentEnd(other.mAlignmentEnd), mMateAlignmentStart(other.mMateAlignmentStart), mMappingQuality(other.mMappingQuality), mInferredInsertSize(other.mInferredInsertSize),
mReferenceName(other.mReferenceName), mMateReferenceName(other.mMateReferenceName), mReadName(other.mReadName), readGroup(other.readGroup){
    if(other.mAttributes != nullptr)
        mAttributes = other.mAttributes;
    else
        mAttributes = nullptr;

    //---why use deep copy here?
    if(other.mReadBases != nullptr){
        mReadBases = std::shared_ptr<uint8_t[]>(new uint8_t [baseLength+1]{0});
        memcpy(mReadBases.get(), other.mReadBases.get(), baseLength);
    }
    else{
        mReadBases = nullptr;
    }
    if(other.mBaseQualities != nullptr) {
        mBaseQualities = std::shared_ptr<uint8_t[]>(new uint8_t[baseQualitiesLength+1]{0});
        memcpy(mBaseQualities.get(), other.mBaseQualities.get(), baseQualitiesLength);
    }
    else{
        mBaseQualities = nullptr;
    }
    mCigar = std::make_shared<Cigar>(other.mCigar->getCigarElements());
}

std::shared_ptr<SimpleInterval> SAMRecord::getLoc() {
    return std::make_shared<SimpleInterval>(ContigMap::getContigInt(mReferenceName), mAlignmentStart, mAlignmentEnd);
}

int SAMRecord::getEndAfterFliter() const {
    return mAlignmentEnd;
}

PositionToCigar::PositionToCigar(int cigarOffset, int currentStart, int offset):cigarOffset(cigarOffset), currentStart(currentStart), offset(offset) {}

bool SAMRecord::overlaps(std::shared_ptr<Locatable> other) {
    return overlaps(getStart(), getEnd(), other->getStart(), other->getEnd());
}

bool SAMRecord::overlaps(int start, int end, int start2, int end2) {
    return (start2 >= start && start2 <= end) || (end2 >=start && end2 <= end) || encloses(start2, end2, start, end);
}

bool SAMRecord::encloses(int outerStart, int outerEnd, int innerStart, int innerEnd) {
    return innerStart >= outerStart && innerEnd <= outerEnd;
}
