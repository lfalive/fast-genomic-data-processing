//
// Created by 梦想家xixi on 2021/11/10.
//

#ifndef MUTECT2CPP_MASTER_ALIGNMENTUTILS_H
#define MUTECT2CPP_MASTER_ALIGNMENTUTILS_H

#include <memory>
#include <Haplotype.h>
#include "cigar/Cigar.h"
#include "samtools/SAMRecord.h"
#include "SmithWatermanAligner.h"

class CigarPairTransform;

class AlignmentUtils {
public:
    static std::vector<CigarPairTransform> cigarPairTransformers;
    static std::set<CigarOperator> ALIGNED_TO_GENOME_PLUS_SOFTCLIPS;
    static std::shared_ptr<Cigar> consolidateCigar(std::shared_ptr<Cigar> c);
    static bool needsConsolidation(const std::shared_ptr<Cigar>& c);

    /**
     * Get the byte[] from bases that cover the reference interval refStart -> refEnd given the
     * alignment of bases to the reference (basesToRefCigar) and the start offset of the bases on the reference
     *
     * refStart and refEnd are 0 based offsets that we want to obtain.  In the client code, if the reference
     * bases start at position X and you want Y -> Z, refStart should be Y - X and refEnd should be Z - X.
     *
     * If refStart or refEnd would start or end the new bases within a deletion, this function will return null
     *
     * @param bases
     * @param refStart
     * @param refEnd
     * @param basesStartOnRef where does the bases array start w.r.t. the reference start?  For example, bases[0] of
     *                        could be at refStart == 0 if basesStartOnRef == 0, but it could just as easily be at
     *                        10 (meaning bases doesn't fully span the reference), which would be indicated by basesStartOnRef == 10.
     *                        It's not trivial to eliminate this parameter because it's tied up with the cigar
     * @param basesToRefCigar the cigar that maps the bases to the reference genome
     * @return a byte[] containing the bases covering this interval, or null if we would start or end within a deletion
     */
    static std::pair<int, std::shared_ptr<uint8_t[]>> getBasesCoveringRefInterval(int refStart, int refEnd, std::shared_ptr<uint8_t[]> bases, int length, int basesStartOnRef, const std::shared_ptr<Cigar>& basesToRefCigar);

    static std::shared_ptr<Cigar> trimCigarByReference(const std::shared_ptr<Cigar>& cigar, int start, int end);

    /**
     * Does cigar start or end with a deletion operation?
     *
     * @param cigar a non-null cigar to test
     * @return true if the first or last operator of cigar is a D
     */
    static bool startsOrEndsWithInsertionOrDeletion(const std::shared_ptr<Cigar>& cigar);

    /**
     * Removing a trailing deletion from the incoming cigar if present
     *
     * @param c the cigar we want to update
     * @return a non-null Cigar
     */
    static std::shared_ptr<Cigar> removeTrailingDeletions(std::shared_ptr<Cigar> c);

    static std::shared_ptr<Cigar> trimCigarByBases(const std::shared_ptr<Cigar>& cigar, int start, int end);

    static std::shared_ptr<Cigar> leftAlignSingleIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex, bool cleanupCigar);

    static std::shared_ptr<Cigar> leftAlignSingleIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex, int leftmostAllowedAlignment, bool cleanupCigar1);

    /**
     * Takes the alignment of the read sequence <code>readSeq</code> to the reference sequence <code>refSeq</code>
     * starting at 0-based position <code>refIndex</code> on the <code>refSeq</code> and specified by its <code>cigar</code>.
     * The last argument <code>readIndex</code> specifies 0-based position on the read where the alignment described by the
     * <code>cigar</code> starts. Usually cigars specify alignments of the whole read to the ref, so that readIndex is normally 0.
     * Use non-zero readIndex only when the alignment cigar represents alignment of a part of the read. The refIndex in this case
     * should be the position where the alignment of that part of the read starts at. In other words, both refIndex and readIndex are
     * always the positions where the cigar starts on the ref and on the read, respectively.
     * <p/>
     * If the alignment has one or more indels, this method attempts to move them left across a stretch of repetitive bases.
     * For instance, if the original cigar specifies that (any) one AT is deleted from a repeat sequence TATATATA, the output
     * cigar will always mark the leftmost AT as deleted. If there is no indel in the original cigar or if the indel position
     * is determined unambiguously (i.e. inserted/deleted sequence is not repeated), the original cigar is returned.
     *
     * Note that currently we do not actually support the case where there is more than one indel in the alignment.  We will throw
     * an exception if there is -- unless the
     *
     * @param cigar     structure of the original alignment
     * @param refSeq    reference sequence the read is aligned to
     * @param readSeq   read sequence
     * @param refIndex  0-based alignment start position on ref
     * @param readIndex 0-based alignment start position on read
     * @param leftmostAllowedAlignment left align indel no further left than this index (0-based)
     * @param doNotThrowExceptionForMultipleIndels  if true we will not throw an exception if we encounter multiple indels in the alignment will instead will return the original cigar
     * @return a non-null cigar, in which the indels are guaranteed to be placed at the leftmost possible position across a repeat (if any)
     */
    static std::shared_ptr<Cigar> leftAlignIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex, int leftmostAllowedAlignment,  bool doNotThrowExceptionForMultipleIndels);

    static std::shared_ptr<Cigar> leftAlignIndel(std::shared_ptr<Cigar> cigar, std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex, bool doNotThrowExceptionForMultipleIndels);

    static std::shared_ptr<Cigar> cleanUpCigar(const std::shared_ptr<Cigar>& c);

    static std::shared_ptr<SAMRecord> createReadAlignedToRef(const std::shared_ptr<SAMRecord>& originalRead, const std::shared_ptr<Haplotype>& haplotype,
                                                             const std::shared_ptr<Haplotype>& refHaplotype, int referenceStart, bool isInformative,
                                                             SmithWatermanAligner *aligner);

    /**
    * Is the offset inside a deletion?
    *
    * @param cigar         the read's CIGAR -- cannot be null
    * @param offset        the offset into the CIGAR
    * @return true if the offset is inside a deletion, false otherwise
    */
    static bool isInsideDeletion(std::shared_ptr<Cigar> cigar, int offset);

    /**
     * Calculate the index into the read's bases of the beginning of the encompassing cigar element for a given cigar and offset
     *
     * @param cigar            the read's CIGAR -- cannot be null
     * @param offset           the offset to use for the calculation or -1 if in the middle of a deletion
     * @param isDeletion       are we in the middle of a deletion?
     * @param alignmentStart   the alignment start of the read
     * @param refLocus         the reference position of the offset
     * @return a non-negative int index
     */
    static int calcAlignmentByteArrayOffset( std::shared_ptr<Cigar> cigar, int offset, bool isDeletion, int alignmentStart, int refLocus);

    /**
    * Get the number of bases aligned to the genome, including soft clips
    *
    * If read is not mapped (i.e., doesn't have a cigar) returns 0
    *
    * @param r a non-null Read
    * @return the number of bases aligned to the genome in R, including soft clipped bases
    */
    static int getNumAlignedBasesCountingSoftClips(std::shared_ptr<SAMRecord> r);

private:
    static std::shared_ptr<Cigar> trimCigar(const std::shared_ptr<Cigar>& cigar, int start, int end, bool byReference);

    static void ensureLeftAlignmentHasGoodArguments(const std::shared_ptr<Cigar>& cigar, std::shared_ptr<uint8_t[]> refSeq, std::shared_ptr<uint8_t[]> readSeq, int refIndex, int readIndex);

    static std::shared_ptr<uint8_t[]> createIndelString(const std::shared_ptr<Cigar>& cigar, int indexOfIndel, std::shared_ptr<uint8_t[]> refSeq, int refLength, std::shared_ptr<uint8_t[]> readSeq, int readLength, int refIndex, int readIndex, int &);

    static std::shared_ptr<Cigar> moveCigarLeft(const std::shared_ptr<Cigar>& cigar, int indexOfIndel);

    /**
     * Counts the number of I/D operators
     *
     * @param cigar   cigar to check -- cannot be null
     * @return  non-negative count of indel operators
     */
    static int countIndelElements( std::shared_ptr<Cigar>& cigar);

protected:
    /**
     * Helper function for trimCigar that adds cigar elements (of total length X) of elt.op to dest for
     * X bases that fall between start and end, where the last position of the base is pos.
     *
     * The primary use of this function is to create a new cigar element list that contains only
     * elements that occur between start and end bases in an initial cigar.
     *
     * Note that this function may return multiple cigar elements (1M1M etc) that are best consolidated
     * after the fact into a single simpler representation.
     *
     * @param dest we will append our cigar elements to this list
     * @param pos the position (0 indexed) where elt started
     * @param start only include bases that occur >= this position
     * @param end only include bases that occur <= this position
     * @param elt the element we are slicing down
     * @return the position after we've traversed all elt.length bases of elt
     */
    static int addCigarElements(std::vector<CigarElement> & dest, int pos, int start, int end, CigarElement elt);

    static bool isIndelAlignedTooFarLeft(const std::shared_ptr<Cigar>& cigar, int leftmostAllowedAlignment);

    static bool cigarHasZeroSizeElement(const std::shared_ptr<Cigar>& c);

    /**
     * Get the offset (base 0) of the first reference aligned base in Cigar that occurs after readStartByBaseOfCigar base of the cigar
     *
     * The main purpose of this routine is to find a good start position for a read given it's cigar.  The real
     * challenge is that the starting base might be inside an insertion, in which case the read actually starts
     * at the next M/EQ/X operator.
     *
     * @param cigar a non-null cigar
     * @param readStartByBaseOfCigar finds the first base after this (0 indexed) that aligns to the reference genome (M, EQ, X)
     * @return an offset into cigar
     */
    static int calcFirstBaseMatchingReferenceInCigar(const std::shared_ptr<Cigar>& cigar, int readStartByBaseOfCigar);

    /**
    * Generate a new Cigar that maps the operations of the first cigar through those in a second
    *
    * For example, if first is 5M and the second is 2M1I2M then the result is 2M1I2M.
    * However, if first is 1M2D3M and second is 2M1I3M this results in a cigar X
    *
    * ref   : AC-GTA
    * hap   : ACxGTA  - 2M1I3M
    * read  : A--GTA  - 1M2D3M
    * result: A--GTA => 1M1D3M
    *
    * ref   : ACxG-TA
    * hap   : AC-G-TA  - 2M1D3M
    * read  : AC-GxTA  - 3M1I2M
    * result: AC-GxTA => 2M1D1M1I2M
    *
    * ref   : ACGTA
    * hap   : ACGTA  - 5M
    * read  : A-GTA  - 1M1I3M
    * result: A-GTA => 1M1I3M
    *
    * ref   : ACGTAC
    * hap   : AC---C  - 2M3D1M
    * read  : AC---C  - 3M
    * result: AG---C => 2M3D
    *
    * The constraint here is that both cigars should imply that the result have the same number of
    * reference bases (i.e.g, cigar.getReferenceLength() are equals).
    *
    * @param firstToSecond the cigar mapping hap1 -> hap2
    * @param secondToThird the cigar mapping hap2 -> hap3
    * @return A cigar mapping hap1 -> hap3
    */
    static std::shared_ptr<Cigar> applyCigarToCigar(const std::shared_ptr<Cigar>& firstToSecond, const std::shared_ptr<Cigar>& secondToThird);

    static CigarPairTransform getTransformer(CigarOperator op12, CigarOperator op23);
};

/**
 * transformations that project one alignment state through another
 *
 * Think about this as a state machine, where we have:
 *
 * bases3 : xxx A zzz
 * bases2 : xxx B zzz
 * bases1 : xxx C zzz
 *
 * where A, B and C are alignment states of a three way alignment.  We want to capture
 * the transition from operation mapping 1 -> 2 and an operation mapping 2 -> 3 and its
 * associated mapping from 1 -> 3 and the advancement of the cigar states of 1->2 and 2->3.
 *
 * Imagine that A, B, and C are all equivalent (so that op12 = M and op23 = M).  This implies
 * a mapping of 1->3 of M, and in this case the next states to consider in the 3 way alignment
 * are the subsequent states in 1 and 2 (so that advance12 and advance23 are both 1).
 *
 * Obviously not all of the states and their associated transitions are so simple.  Suppose instead
 * that op12 = I, and op23 = M.  What does this look like:
 *
 * bases3 : xxx - A zzz
 * bases2 : xxx - B zzz
 * bases1 : xxx I C zzz
 *
 * It means that op13 must be an insertion (as we have an extra base in 1 thats not present in 2 and
 * so not present in 3).  We advance the cigar in 1 by 1 (as we've consumed one base in 1 for the I)
 * but we haven't yet found the base corresponding to the M of op23.  So we don't advance23.
 */
class CigarPairTransform {
private:
    static std::set<enum CigarOperator> getCigarSet(CigarOperator masterOp){
        switch (masterOp) {
            case M:
                return std::set<enum CigarOperator>{CigarOperator::M, CigarOperator::EQ, CigarOperator::X};
            case I:
                return std::set<enum CigarOperator>{CigarOperator::I, CigarOperator::S};
            case D:
                return std::set<enum CigarOperator>{CigarOperator::D};
            default:
                throw std::exception();
        }
    }

public:
    std::set<enum CigarOperator> op12, op23;
    CigarOperator op13;
    int advance12, advance23;

    CigarPairTransform(CigarOperator op12, CigarOperator op23, CigarOperator op13, int advance12, int advance23 ){
        this->op12 = getCigarSet(op12);
        this->op23 = getCigarSet(op23);
        this->op13 = op13;
        this->advance12 = advance12;
        this->advance23 = advance23;
    }
};

#endif //MUTECT2CPP_MASTER_ALIGNMENTUTILS_H
