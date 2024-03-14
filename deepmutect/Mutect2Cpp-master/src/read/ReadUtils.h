//
// Created by 梦想家xixi on 2021/12/18.
//

#ifndef MUTECT2CPP_MASTER_READUTILS_H
#define MUTECT2CPP_MASTER_READUTILS_H

#include "samtools/SAMRecord.h"
#include "samtools/SAMFileHeader.h"

enum ClippingTail {
    LEFT_TAIL,
    RIGHT_TAIL
};

enum Passes {
    FIRST,
    SECOND,
    END
};



class ReadUtils {
public:
    static const int CANNOT_COMPUTE_ADAPTOR_BOUNDARY = INT32_MIN;
    static const int CLIPPING_GOAL_NOT_REACHED;

    static std::shared_ptr<SAMRecord> emptyRead(SAMRecord & read);
    static std::shared_ptr<SAMRecord> emptyRead(std::shared_ptr<SAMRecord> & read);
    static void assertAttributeNameIsLegal(std::string& attributeName);
    static int getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> & read, int refCoord, ClippingTail tail);
    static int getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord> & read, int refCoord, ClippingTail tail, bool allowGoalNotReached);
    static int getReadCoordinateForReferenceCoordinate(int alignmentStart, std::shared_ptr<Cigar> cigar, int refCoord, ClippingTail tail, bool allowGoalNotReached);
    static CigarElement* readStartsWithInsertion(const std::shared_ptr<Cigar>& cigarForRead);
    static CigarElement* readStartsWithInsertion(const std::shared_ptr<Cigar>& cigarForRead, bool ignoreSoftClipOps);
    static int getSoftStart(std::shared_ptr<SAMRecord> & read);
    static int getSoftEnd(std::shared_ptr<SAMRecord> & read);
    static bool hasBaseIndelQualities(std::shared_ptr<SAMRecord> & read);
    const static char DEFAULT_INSERTION_DELETION_QUAL = 45;
    static std::string BQSR_BASE_INSERTION_QUALITIES;
    static std::string BQSR_BASE_DELETION_QUALITIES;
    static std::shared_ptr<uint8_t[]> getExistingBaseInsertionQualities(SAMRecord& read, int & length);
    static std::shared_ptr<uint8_t[]> getExistingBaseDeletionQualities(SAMRecord& read, int & length);
    static std::shared_ptr<uint8_t[]> getExistingBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static std::shared_ptr<uint8_t[]> getExistingBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static std::shared_ptr<uint8_t[]> getBaseInsertionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static std::shared_ptr<uint8_t[]> getBaseDeletionQualities(std::shared_ptr<SAMRecord> & read, int & length);
    static std::shared_ptr<uint8_t[]> getBaseInsertionQualities(SAMRecord & read, int & length);
    static std::shared_ptr<uint8_t[]> getBaseDeletionQualities(SAMRecord & read, int & length);
    static void setInsertionBaseQualities(std::shared_ptr<SAMRecord> & read, std::shared_ptr<uint8_t[]> quals, int length);
    static void setDeletionBaseQualities(std::shared_ptr<SAMRecord> & read, std::shared_ptr<uint8_t[]> quals, int length);
    static int getAssignedReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader* header);
    static int getSAMFlagsForRead(std::shared_ptr<SAMRecord> & read);
    static int getMateReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader* header);
    static bool alignmentAgreesWithHeader(SAMFileHeader* header, std::shared_ptr<SAMRecord> & read);
    static int getReferenceIndex(std::shared_ptr<SAMRecord> & read, SAMFileHeader* header);
    static bool hasWellDefinedFragmentSize(std::shared_ptr<SAMRecord> & read);
    static bool hasWellDefinedFragmentSize(SAMRecord * read);
    static int getAdaptorBoundary(SAMRecord * read);
    static bool isBaseInsideAdaptor(std::shared_ptr<SAMRecord> & read, long basePos);
    static bool isInsideRead(std::shared_ptr<SAMRecord> & read, int referenceCoordinate);
    static bool isF2R1(std::shared_ptr<SAMRecord> read);

    // overrite some method without using SAMRecord class for performance
    static bool hasWellDefinedFragmentSize(bam1_t * read, sam_hdr_t * hdr);
    static int getAdaptorBoundary(bam1_t * read, sam_hdr_t * hdr);
    static const char * getReferenceName(bam1_t * read, sam_hdr_t * hdr);
    static bool isPaired(bam1_t * read);
    static bool isProperlyPaired(bam1_t * read);
    static bool isUnmapped(bam1_t * read, sam_hdr_t * hdr);
    static bool isReverseStrand(bam1_t * read);
    static bool getMateUnmappedFlagUnchecked(bam1_t * read);
    static bool getMateNegativeStrandFlagUnchecked(bam1_t * read);
    static bool getProperPairFlagUnchecked(bam1_t * read);
    static bool mateIsUnmapped(bam1_t * read, sam_hdr_t * hdr);
    static bool mateIsReverseStrand(bam1_t * read);
    static bool alignmentAgreesWithHeader(sam_hdr_t * hdr, bam1_t * read);
    static hts_pos_t getEnd(bam1_t * read);

    static uint8_t decodeBase(uint8_t base);

    static bool consumesReadBases(uint32_t cigarElement);
    static bool consumesReferenceBases(uint32_t cigarElement);


    /**
     * Calculates how much the alignment should be shifted when hard/soft clipping is applied
     * to the cigar
     */
    static int calculateAlignmentStartShift(int n_cigar, uint32_t * oldCigar, int newReadBasesClipped);

    /**
     * Returns the read coordinate corresponding to the requested reference coordinate.
     *
     * WARNING: if the requested reference coordinate happens to fall inside or just before a deletion (or skipped region) in the read, this function
     * will return the last read base before the deletion (or skipped region). This function returns a
     * Pair(int readCoord, boolean fallsInsideOrJustBeforeDeletionOrSkippedRegion) so you can choose which readCoordinate to use when faced with
     * a deletion (or skipped region).
     *
     * SUGGESTION: Use getReadCoordinateForReferenceCoordinate(GATKSAMRecord, int, ClippingTail) instead to get a
     * pre-processed result according to normal clipping needs. Or you can use this function and tailor the
     * behavior to your needs.
     *
     * @param read
     * @param refCoord the requested reference coordinate
     * @return the read coordinate corresponding to the requested reference coordinate. (see warning!)
     */
    static std::pair<int, bool> getReadCoordinateForReferenceCoordinate(std::shared_ptr<SAMRecord>& read, int refCoord);

private:
    static std::pair<int, bool> getReadCoordinateForReferenceCoordinate(int alignmentStart, std::shared_ptr<Cigar> cigar, int refCoord, bool allowGoalNotReached);

};


#endif //MUTECT2CPP_MASTER_READUTILS_H
