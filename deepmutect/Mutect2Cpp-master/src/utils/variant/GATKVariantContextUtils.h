//
// Created by lhh on 4/26/22.
//

#ifndef MUTECT2CPP_MASTER_GATKVARIANTCONTEXTUTILS_H
#define MUTECT2CPP_MASTER_GATKVARIANTCONTEXTUTILS_H


#include <cstdint>
#include <VariantContext.h>
#include "Locatable.h"
#include "parallel_hashmap/phmap.h"

enum FilteredRecordMergeType {
    /**
     * Union - leaves the record if any record is unfiltered.
     */
    KEEP_IF_ANY_UNFILTERED,
    /**
     * Requires all records present at site to be unfiltered. VCF files that don't contain the record don't influence this.
     */
    KEEP_IF_ALL_UNFILTERED,
    /**
     * If any record is present at this site (regardless of possibly being filtered), then all such records are kept and the filters are reset.
     */
    KEEP_UNCONDITIONAL
};

enum GenotypeMergeType {
    /**
     * Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
     */
    UNIQUIFY,
    /**
     * Take genotypes in priority order (see the priority argument).
     */
    PRIORITIZE,
    /**
     * Take the genotypes in any order.
     */
    UNSORTED,
    /**
     * Require that all samples/genotypes be unique between all inputs.
     */
    REQUIRE_UNIQUE
};

class AlleleMapper {
private:
    std::shared_ptr<VariantContext> vc;
    std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>> map;
public:
    AlleleMapper(std::shared_ptr<VariantContext> vc);
    AlleleMapper(std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>> map);

    bool needsRemapping();

    std::shared_ptr<Allele> remap(std::shared_ptr<Allele>& a);
    std::shared_ptr<std::vector<std::shared_ptr<Allele>>> remap(std::vector<std::shared_ptr<Allele>> & as);
    std::shared_ptr<std::vector<std::shared_ptr<Allele>>> values();
};

class GATKVariantContextUtils {
private:
    static bool compareVariantContext(std::shared_ptr<VariantContext>& vc1, std::shared_ptr<VariantContext>& vc2, phmap::flat_hash_map<std::string, int>& ComparatorMap);


    static std::shared_ptr<Allele> determineReferenceAllele(std::vector<std::shared_ptr<VariantContext>> & VCs);

    /**
    * Determines the common reference allele
    *
    * @param VCs    the list of VariantContexts
    * @param loc    if not null, ignore records that do not begin at this start location
    * @return possibly null Allele
    */
    static std::shared_ptr<Allele> determineReferenceAllele(std::vector<std::shared_ptr<VariantContext>> & VCs, Locatable* loc);

    static bool isNonSymbolicExtendableAllele(std::shared_ptr<Allele>& allele);

    static void mergeGenotypes(GenoTypesContext& mergedGenotypes, std::shared_ptr<VariantContext> oneVC, std::shared_ptr<AlleleMapper> alleleMapping, bool uniquifySamples);

    static std::string mergedSampleName(std::string trackName, std::string  sampleName, bool uniquify);

    static bool hasPLIncompatibleAlleles(std::vector<std::shared_ptr<Allele>>& alleleSet1, std::vector<std::shared_ptr<Allele>> & alleleSet2);

public:
    /**
     * Finds number of repetitions a string consists of.
     * Same as {@link #findNumberOfRepetitions} but operates on subarrays of a bigger array to save on copying.
     * For example, for string ATAT and repeat unit AT, number of repetitions = 2
     * @param repeatUnitFull             Non-empty substring represented by byte array
     * @param repeatUnitFullLength       the total length of repeatUnitFull
     * @param offsetInRepeatUnitFull     the offset in repeatUnitFull from which to read the repeat unit
     * @param repeatUnitLength           length of the repeat unit
     * @param testStringFull             string to test (represented by byte array), may be empty
     * @param testStringFullLength       the total length of testStringFull
     * @param offsetInTestStringFull     the offset in offsetInRepeatUnitFull from which to read the test string
     * @param testStringLength           length of the test string
     * @param leadingRepeats         Look for leading (at the beginning of string) or trailing (at end of string) repetitions
     * For example:
     *    GATAT has 0 leading repeats of AT but 2 trailing repeats of AT
     *    ATATG has 1 leading repeat of A but 2 leading repeats of AT
     *    CCCCCCCC has 2 leading and 2 trailing repeats of CCC
     * @return  Number of repetitions (0 if testString is not a concatenation of n repeatUnit's, including the case of empty testString)
     */
    static int findNumberOfRepetitions(uint8_t* repeatUnitFull, int repeatUnitFullLength, int offsetInRepeatUnitFull, int repeatUnitLength, uint8_t* testStringFull, int testStringFullLength, int offsetInTestStringFull, int testStringLength, bool leadingRepeats);

    static int findNumberOfRepetitions(uint8_t* repeatUnit, int repeatUnitLength, uint8_t* testString, int testStringLength, bool leadingRepeats);

    /**
    * Merges VariantContexts into a single hybrid.  Takes genotypes for common samples in priority order, if provided.
    * If uniquifySamples is true, the priority order is ignored and names are created by concatenating the VC name with
    * the sample name
    *
    * @param unsortedVCs               collection of unsorted VCs
    * @param priorityListOfVCs         priority list detailing the order in which we should grab the VCs
    * @param filteredRecordMergeType   merge type for filtered records
    * @param genotypeMergeOptions      merge option for genotypes
    * @param filteredAreUncalled       are filtered records uncalled?
    * @return new VariantContext       representing the merge of unsortedVCs
    */
    static std::shared_ptr<VariantContext> simpleMerge(std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> unsortedVCs, std::vector<std::string>& priorityListOfVCs,  FilteredRecordMergeType filteredRecordMergeType,
                                                GenotypeMergeType genotypeMergeOptions, bool filteredAreUncalled);

    static std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> sortVariantContextsByPriority(std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> unsortedVCs, std::vector<std::string>& priorityListOfVCs, GenotypeMergeType genotypeMergeOptions);

    static bool contextMatchesLoc(std::shared_ptr<VariantContext>& vc, Locatable* loc);

    static std::shared_ptr<AlleleMapper> resolveIncompatibleAlleles(std::shared_ptr<Allele> refAllele, std::shared_ptr<VariantContext> vc, std::vector<std::shared_ptr<Allele>>& alleles);

    /**
    * Create an allele mapping for the given context where its reference allele must (potentially) be extended to the given allele
    *
    * The refAllele is the longest reference allele seen at this start site.
    * So imagine it is:
    * refAllele: ACGTGA
    * myRef:     ACGT
    * myAlt:     A
    *
    * We need to remap all of the alleles in vc to include the extra GA so that
    * myRef => refAllele and myAlt => AGA
    *
    * @param refAllele          the new (extended) reference allele
    * @param oneVC              the Variant Context to extend
    * @param currentAlleles     the list of alleles already created
    * @return a non-null mapping of original alleles to new (extended) ones
    */
    static std::shared_ptr<std::map<std::shared_ptr<Allele>, std::shared_ptr<Allele>>> createAlleleMapping(std::shared_ptr<Allele> refAllele, std::shared_ptr<VariantContext> oneVc, phmap::flat_hash_set<std::shared_ptr<Allele>> &currentAlleles);

    static std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>> createAlleleMapping(std::shared_ptr<Allele> refAllele, std::shared_ptr<VariantContext> oneVc, const std::vector<std::shared_ptr<Allele>> &currentAlleles);

    static std::shared_ptr<GenoTypesContext> stripPLsAndAD(std::shared_ptr<GenoTypesContext> genotypes);

    static std::shared_ptr<Genotype> removePLsAndAD(std::shared_ptr<Genotype> g);

    /**
    * Trim the alleles in inputVC forward and reverse, as requested
    *
    * @param inputVC a non-null input VC whose alleles might need a haircut
    * @param trimForward should we trim up the alleles from the forward direction?
    * @param trimReverse should we trim up the alleles from the reverse direction?
    * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
    */
    static std::shared_ptr<VariantContext> trimAlleles(std::shared_ptr<VariantContext> inputVC,  bool trimForward, bool trimReverse);

    static int computeReverseClipping(std::vector<std::shared_ptr<Allele>>& unclippedAlleles, std::shared_ptr<uint8_t[]> ref, int refLength);

    /**
     * Clip out any unnecessary bases off the front of the alleles
     *
     * The VCF spec represents alleles as block substitutions, replacing AC with A for a
     * 1 bp deletion of the C.  However, it's possible that we'd end up with alleles that
     * contain extra bases on the left, such as GAC/GA to represent the same 1 bp deletion.
     * This routine finds an offset among all alleles that can be safely trimmed
     * off the left of each allele and still represent the same block substitution.
     *
     * A/C => A/C
     * AC/A => AC/A
     * ACC/AC => CC/C
     * AGT/CAT => AGT/CAT
     * <DEL>/C => <DEL>/C
     *
     * @param unclippedAlleles a non-null list of alleles that we want to clip
     * @return the offset into the alleles where we can safely clip, inclusive, or
     *   -1 if no clipping is tolerated.  So, if the result is 0, then we can remove
     *   the first base of every allele.  If the result is 1, we can remove the
     *   second base.
     */
    static int computeForwardClipping(std::vector<std::shared_ptr<Allele>>& unclippedAlleles);

    /**
    * Trim up alleles in inputVC, cutting out all bases up to fwdTrimEnd inclusive and
    * the last revTrim bases from the end
    *
    * @param inputVC a non-null input VC
    * @param fwdTrimEnd bases up to this index (can be -1) will be removed from the start of all alleles
    * @param revTrim the last revTrim bases of each allele will be clipped off as well
    * @return a non-null VariantContext (may be == to inputVC) with trimmed up alleles
    */
    static std::shared_ptr<VariantContext> trimAlleles(std::shared_ptr<VariantContext> inputVC, int fwdTrimEnd, int revTrim);

    static std::shared_ptr<GenoTypesContext> updateGenotypesWithMappedAlleles(std::shared_ptr<GenoTypesContext> originalGenotypes, AlleleMapper& alleleMapper);

    static std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>> getNumTandemRepeatUnits(const std::shared_ptr<VariantContext>& vc, std::shared_ptr<uint8_t[]> refBasesStartingAtVCWithPad, int len);

    static std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>> getNumTandemRepeatUnits(const std::shared_ptr<uint8_t[]> & refBases, int refLen, const std::shared_ptr<uint8_t[]> & altBases, int altLen,
                                                                                     const std::shared_ptr<uint8_t[]> & remainingRefContext, int remainLen);

    static int findRepeatedSubstring(uint8_t * bases, int basesLen);
};


#endif //MUTECT2CPP_MASTER_GATKVARIANTCONTEXTUTILS_H
