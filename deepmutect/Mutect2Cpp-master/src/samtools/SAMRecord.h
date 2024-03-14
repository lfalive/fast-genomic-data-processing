//
// Created by 梦想家xixi on 2021/12/20.
//

#ifndef MUTECT2CPP_MASTER_SAMRECORD_H
#define MUTECT2CPP_MASTER_SAMRECORD_H

#include <string>
#include "Locatable.h"
#include "Cigar.h"
#include "SAMBinaryTagAndValue.h"
#include "htslib/sam.h"
#include "SAMFileHeader.h"
#include "SimpleInterval.h"

struct PositionToCigar{
    PositionToCigar(int cigarOffset, int currentStart, int offset);

    int cigarOffset;
    int currentStart;
    int offset;
};

class SAMRecord{
private:
    std::shared_ptr<uint8_t[]> mReadBases;
    int baseLength;
    uint8_t readGroup;
    std::shared_ptr<uint8_t[]> mBaseQualities;
    int baseQualitiesLength;
    std::string mReadName;
    std::string mReferenceName;
    int mAlignmentStart;
    int mAlignmentEnd;
    int mMappingQuality;
    int adaptorBoundary = -1;
    bool isCalAdaptorBoundary = false;
    std::shared_ptr<Cigar> mCigar;
    int mFlags;
    std::string mMateReferenceName;
    int mMateAlignmentStart;
    int mInferredInsertSize;
    std::shared_ptr< SAMBinaryTagAndValue> mAttributes;



    void setFlag(bool flag, int bit);
    void requireReadPaired() const;
    bool getMateUnmappedFlagUnchecked() const;

    /**
    * Checks to see if the two sets of coordinates have any overlap.
    */
    static bool overlaps(int start, int end, int start2, int end2);

    static bool encloses(int outerStart, int outerEnd, int innerStart, int innerEnd);

public:
    // map from position to cigar index
    std::vector<std::pair<int, PositionToCigar>> PositionToCigarMap;


    SAMRecord(std::shared_ptr<uint8_t[]> base, int baseLength, std::shared_ptr<uint8_t[]> baseQualities, int baseQualitiesLength, std::string &name);
    SAMRecord(bam1_t * read, sam_hdr_t * hdr, bool load = true);
    SAMRecord(const SAMRecord & other);
    ~SAMRecord();
    static const std::string NO_ALIGNMENT_REFERENCE_NAME;
    static const int NO_ALIGNMENT_START = 0;
    static const int NO_MAPPING_QUALITY = 0;
    static const int NO_ALIGNMENT_REFERENCE_INDEX = -1;
    std::string& getName();
    void setName(std::string& name);
    int getLength() const;
    bool isEmpty() const {return getLength() == 0;}
    void setPosition(std::string& contig, int start);
    void setPosition(Locatable* locatable);
    std::string & getAssignedContig();
    int getAssignedStart() const;
    int getUnclippedStart();
    int getUnclippedEnd();
    bool getReadUnmappedFlag() const;
    bool isUnmapped();
    bool mateIsUnmapped();
    bool getMateUnmappedFlag();
    bool isPaired() const;
    std::string& getMateContig();
    int getMateStart();
    void setMatePosition(std::string& contig, int start);
    void setIsPaired(bool isPaired);
    void setReadPairedFlag(bool flag);
    void setProperPairFlag(bool flag);
    void setMatePosition(Locatable* locatable);
    int getFragmentLength() const;
    void setFragmentLength(int fragmentLength);
    int getMappingQuality() const;
    void setMappingQuality(int mappingQuality);
    std::shared_ptr<uint8_t[]> getBases();
    std::shared_ptr<uint8_t[]> & getBasesNoCopy();
    uint8_t getBase(const int i) {return mReadBases[i];}
    void setBases(std::shared_ptr<uint8_t[]>bases, int length);
    std::shared_ptr<uint8_t[]> getBaseQualities();
    std::shared_ptr<uint8_t[]> & getBaseQualitiesNoCopy();
    int getBaseQualitiesLength() const;
    uint8_t getBaseQuality(const int i) {return mBaseQualities[i];}
    void setBaseQualities(std::shared_ptr<uint8_t[]> baseQualities, int length);
    const std::shared_ptr<Cigar> & getCigar();
    std::vector<CigarElement>& getCigarElements();
    CigarElement getCigarElement(int index);
    void setCigar(std::shared_ptr<Cigar> cigar);
    bool getProperPairFlag();
    bool getProperPairFlagUnchecked() const;
    bool isProperlyPaired();
    void setIsProperlyPaired(bool isProperlyPaired);
    int getStart() const;
    int getEnd();
    void setIsUnmapped();
    void setReadUnmappedFlag(bool flag);
    void clearAttributes();
    void* getAttribute(short tag);
    std::string& getReadGroup();
    void setAttribute(std::string& attributeName, const std::string& attributeValue);
    void setAttribute(std::string &tag, void* value, Void_Type type, int length);
    void setAttribute(short tag, void* value, Void_Type type, int length);
    int getSoftStart();
    int getSoftEnd();
    std::string & getContig();
	int getContigInt();
    std::string getAttributeAsString(std::string & attributeName);
    bool isReverseStrand() const;
    bool mateIsReverseStrand();
    bool getMateNegativeStrandFlagUnchecked() const;
    bool isFirstOfPair();
    bool getFirstOfPairFlag();
    bool isSecondOfPair();
    bool getSecondOfPairFlag();
    bool isSecondaryAlignment() const;
    bool failsVendorQualityCheck() const;
    bool isDuplicate() const;
    bool isSupplementaryAlignment() const;
    int getAdaptorBoundary();
    int getEndAfterFliter() const;
    std::shared_ptr<SimpleInterval> getLoc();
    void setGroup(uint8_t i) {readGroup = i;}
    uint8_t getGroup() const {return readGroup;}
    bool overlaps(std::shared_ptr<Locatable> other);


private:
    void setAttribute(short tag, void* value, Void_Type type, int length, bool isUnsignedArray);
};


#endif //MUTECT2CPP_MASTER_SAMRECORD_H
