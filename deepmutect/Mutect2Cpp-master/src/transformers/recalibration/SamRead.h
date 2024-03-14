/**
  * A wrapper for sam record
  */

#ifndef SAM_READ_H
#define SAM_READ_H


#include "htslib/sam.h"

//#define PRESERVE_QSCORES_LESS_THAN 6

typedef struct BaseArgument{
    int LowQualityTail;
    int MaximumCycleValue;
    int MismatchesContextSize;
    int MismatchesDefaultQuality;
    int QuantizingLevels;
} BaseArgument;

/**
 * Only to encapsulate sam_hdr_tid2name function in sam.h
 */
const char* GetReferenceName(const sam_hdr_t *h, bam1_t *read);

/**
 * zero-based start
 */
hts_pos_t GetStart(bam1_t *read);

/**
 * @return zero-based end
 */
int GetEnd(bam1_t *read);

int getMateStart(const sam_hdr_t *h, bam1_t *read);

/**
 * Does the read have a position assigned to it for sorting purposes.
 * @return `true if this read has no assigned position or contig.
 */
bool IsUnmapped(const sam_hdr_t *h, bam1_t *read);

/**
 * @return True if this read's mate is unmapped (this includes mates that have a position but are explicitly marked as unmapped,
 *         as well as mates that lack a fully-defined position but are not explicitly marked as unmapped). Otherwise false.
 * @throws IllegalStateException if the read is not paired (has no mate)
 */
bool mateIsUnmapped(const sam_hdr_t *h, bam1_t *read);

/**
 * Returns the observed length of the read's fragment (equivalent to TLEN in SAM).
 *
 * Warning: the precise meaning of this field is implementation/technology dependent.
 *
 * @return The observed length of the fragment (equivalent to TLEN in SAM), or 0 if unknown.
 *         Negative if the mate maps to a lower position than the read.
 */
int getFragmentLength(bam1_t * read);

bool isPaired(bam1_t * read);

/**
 * Can the adaptor sequence of read be reliably removed from the read based on the alignment of
 * read and its mate?
 *
 * @param read the read to check
 * @return true if it can, false otherwise
 */
bool hasWellDefinedFragmentSize(const sam_hdr_t *h, bam1_t * read);

bool getReadPairedFlag(bam1_t * read);

bool getProperPairFlagUnchecked(bam1_t * read);

bool getMateUnmappedFlag(bam1_t *read);

bool getMateUnmappedFlagUnchecked(bam1_t *read);

const char* getMateReferenceName(const sam_hdr_t *h, bam1_t *read);

int getMateAlignmentStart(bam1_t * read);

/**
 * @return True if this read is on the reverse strand as opposed to the forward strand, otherwise false.
 */
bool isReverseStrand(bam1_t *read);

bool mateIsReverseStrand(bam1_t *read);

/**
 * Is a base inside a read?
 *
 * @param read                the read to evaluate
 * @param referenceCoordinate the reference coordinate of the base to test
 * @return true if it is inside the read, false otherwise.
 */
bool isInsideRead(bam1_t *read, int referenceCoordinate);

uint8_t * getBaseQualities(bam1_t * read);

// this function is too expensive
uint8_t getBaseQuality(bam1_t * read, int i);

bool getSecondOfPairFlagUnchecked(bam1_t * read);

/**
 * Calculates the reference coordinate for the beginning of the read taking into account soft clips but not hard clips.
 *
 * Note: getUnclippedStart() adds soft and hard clips, this function only adds soft clips.
 *
 * @return the unclipped start of the read taking soft clips (but not hard clips) into account
 */
int getSoftStart(bam1_t * read);

/**
 * Calculates the reference coordinate for the end of the read taking into account soft clips but not hard clips.
 *
 * Note: getUnclippedEnd() adds soft and hard clips, this function only adds soft clips.
 *
 * @return the unclipped end of the read taking soft clips (but not hard clips) into account
 */
int getSoftEnd(bam1_t * read);

bool consumesReadBases(uint32_t cigarElement);

bool consumesReferenceBases(uint32_t cigarElement);

/**
 * print a human readable format for bases of a read
 */
void PrintBases(bam1_t * read);

void PrintQualities(bam1_t * read);

void PrintCigar(bam1_t * read);

void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar, hts_pos_t *rlen, hts_pos_t *qlen);

#endif