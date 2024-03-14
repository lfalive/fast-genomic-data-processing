//
// Created by 梦想家xixi on 2021/12/8.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H
#define MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H

#include "haplotype/Haplotype.h"
#include "AssemblyRegion.h"
#include "ReferenceCache.h"
#include "Mutect2/AssemblyResultSet.h"
#include "M2ArgumentCollection.h"
#include "ReadThreadingAssembler.h"
#include "PairHMMLikelihoodCalculationEngine.h"
#include "LikelihoodEngineArgumentCollection.h"
#include "SmithWatermanAligner.h"


class AssemblyBasedCallerUtils {
private:
    static std::string phase01;
    static std::string phase10;


    static bool isBiallelic(shared_ptr<VariantContext> vc);

public:
    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    static std::shared_ptr<Haplotype> createReferenceHaplotype(const std::shared_ptr<AssemblyRegion> & region, const std::shared_ptr<SimpleInterval> &referencePadding, ReferenceCache & cache);

    static std::shared_ptr<AssemblyResultSet> assembleReads(const std::shared_ptr<AssemblyRegion>& region, M2ArgumentCollection & argumentCollection, SAMFileHeader* header, ReferenceCache & cache, ReadThreadingAssembler& assemblyEngine);

    static const int REFERENCE_PADDING_FOR_ASSEMBLY = 500;

    static const int MINIMUM_READ_LENGTH_AFTER_TRIMMING = 10;

    static std::shared_ptr<SimpleInterval> getPaddedReferenceLoc(const std::shared_ptr<AssemblyRegion>& region, int referencePadding, SAMFileHeader* header);

    static void finalizeRegion(const std::shared_ptr<AssemblyRegion>& region, bool errorCorrectReads, bool dontUseSoftClippedBases, uint8_t minTailQuality, SAMFileHeader* header, bool correctOverlappingBaseQualities);

    static std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> splitReadsBySample(const std::vector<std::string>& sampleList, const std::string& normalSample, const std::vector<std::shared_ptr<SAMRecord>> & reads);

    /**
     * Instantiates the appropriate likelihood calculation engine.
     */
    static PairHMMLikelihoodCalculationEngine * createLikelihoodCalculationEngine(LikelihoodEngineArgumentCollection& likelihoodArgs);

    /**
     *  Modify base qualities when paired reads overlap to account for the possibility of PCR error.
     *
     *  Overlapping mates provded independent evidence as far as sequencing error is concerned, but their PCR errors
     *  are correlated.  The base qualities are thus limited by the sequencing base quality as well as half of the PCR
     *  quality.  We use half of the PCR quality because downstream we treat read pairs as independent, and summing two halves
     *  effectively gives the PCR quality of the pairs when taken together.
     *
     * @param reads the list of reads to consider
     * @param samplesList   list of samples|
     * @param readsHeader   bam header of reads' source
     * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
     * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR
     * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR
     */
    static void cleanOverlappingReadPairs(vector<shared_ptr<SAMRecord>>& reads, const vector<string>& sampleList, const string& sample, bool setConflictingToZero, int halfOfPcrSnvQual = 0, int halfOfPcrIndelQual = 0);

    // create the assembly using just high quality reads (eg Q20 or higher).  We may want to use lower
    // quality reads in the PairHMM downstream, so we can't use a ReadFilter
    static std::shared_ptr<AssemblyRegion> assemblyRegionWithWellMappedReads(const std::shared_ptr<AssemblyRegion>& originalAssemblyRegion, int minMappingQuality, SAMFileHeader * header);

    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    static shared_ptr<phmap::flat_hash_map<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>> realignReadsToTheirBestHaplotype(AlleleLikelihoods<SAMRecord, Haplotype>& originalReadLikelihoods, shared_ptr<Haplotype>& refHaplotype, shared_ptr<SimpleInterval>& paddedReferenceLoc, SmithWatermanAligner* aligner);

    static double HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY(shared_ptr<Haplotype> h);

    /**
    * Returns the list of events discovered in assembled haplotypes that are active at this location. The results will
    * include events that span the current location if includeSpanningEvents is set to true; otherwise it will only
    * include events that have loc as their start position.
    * @param loc The start position we are genotyping
    * @param haplotypes list of active haplotypes at the current location
    * @param includeSpanningEvents If true, will also return events that span loc
    */
    static shared_ptr<vector<shared_ptr<VariantContext>>> getVariantContextsFromActiveHaplotypes(int loc, vector<shared_ptr<Haplotype>>& haplotypes, bool includeSpanningEvents);

    static shared_ptr<VariantContext> makeMergedVariantContext(shared_ptr<vector<shared_ptr<VariantContext>>> vcs);

    /**
     * Returns a mapping from Allele in the mergedVC, which represents all of the alleles being genotyped at loc,
     * to a list of Haplotypes that support that allele. If the mergedVC includes a spanning deletion allele, all haplotypes
     * that support spanning deletions will be assigned to that allele in the map.
     *
     * @param mergedVC The merged variant context for the locus, which includes all active alternate alleles merged to a single reference allele
     * @param loc The active locus being genotyped
     * @param haplotypes Haplotypes for the current active region
     * @return
     */
    static shared_ptr<phmap::flat_hash_map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Haplotype>>>, hash_Allele, equal_Allele>> createAlleleMapper(shared_ptr<VariantContext> mergedVC, int loc, vector<shared_ptr<Haplotype>>& haplotypes);

    /**
     * Tries to phase the individual alleles based on pairwise comparisons to the other alleles based on all called haplotypes
     *
     * @param calls             the list of called alleles
     * @param calledHaplotypes  the set of haplotypes used for calling
     * @return a non-null list which represents the possibly phased version of the calls
     */
    static shared_ptr<vector<shared_ptr<VariantContext>>> phaseCalls(vector<shared_ptr<VariantContext>>& calls, phmap::flat_hash_set<shared_ptr<Haplotype>>& calledHaplotypes);

    static shared_ptr<map<VariantContext*, shared_ptr<phmap::flat_hash_set<Haplotype*>>>> constructHaplotypeMapping(vector<shared_ptr<VariantContext>>& originalCalls, phmap::flat_hash_set<shared_ptr<Haplotype>>& calledHaplotypes);

    /**
    * Construct the mapping from call (variant context) to phase set ID
    *
    * @param originalCalls    the original unphased calls
    * @param haplotypeMap     mapping from alternate allele to the set of haplotypes that contain that allele
    * @param totalAvailableHaplotypes the total number of possible haplotypes used in calling
    * @param phaseSetMapping  the map to populate in this method;
    *                         note that it is okay for this method NOT to populate the phaseSetMapping at all (e.g. in an impossible-to-phase situation)
    * @return the next incremental unique index
    */
    static int constructPhaseSetMapping(vector<shared_ptr<VariantContext>> &originalCalls, map<VariantContext*, shared_ptr<phmap::flat_hash_set<Haplotype*>>>& haplotypeMap, int totalAvailableHaplotypes, map<VariantContext*, pair<int, string>>& phaseSetMapping);

    static shared_ptr<vector<shared_ptr<VariantContext>>> constructPhaseGroups(vector<shared_ptr<VariantContext>> &originalCalls, map<VariantContext*, pair<int, string>>& phaseSetMapping, int indexTo);

    /**
     * Create a unique identifier given the variant context
     *
     * @param vc   the variant context
     * @return non-null String
     */
    static string createUniqueID(shared_ptr<VariantContext> vc);

    /**
     * Add physical phase information to the provided variant context
     *
     * @param vc   the variant context
     * @param ID   the ID to use
     * @param phaseGT the phase GT string to use
     * @return phased non-null variant context
     */
    static shared_ptr<VariantContext> phaseVC(shared_ptr<VariantContext> vc, string& ID, string& phaseGT, int phaseSetID);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYBASEDCALLERUTILS_H
