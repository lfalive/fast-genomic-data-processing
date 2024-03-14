//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_PAIRHMMLIKELIHOODCALCULATIONENGINE_H
#define MUTECT2CPP_MASTER_PAIRHMMLIKELIHOODCALCULATIONENGINE_H

#include "utils/pairhmm/PairHMM.h"
#include "AssemblyResultSet.h"
#include "utils/genotyper/AlleleLikelihoods.h"
#include "PairHMMNativeArgumentCollection.h"

enum PCRErrorModel{
    /** no specialized PCR error model will be applied; if base insertion/deletion qualities are present they will be used */
    //NONE = 0,
    /** a most aggressive model will be applied that sacrifices true positives in order to remove more false positives */
    HOSTILE = 1,
    /** a more aggressive model will be applied that sacrifices true positives in order to remove more false positives */
    AGGRESSIVE = 2,
    /** a less aggressive model will be applied that tries to maintain a high true positive rate at the expense of allowing more false positives */
    CONSERVATIVE = 3
};

/*
 * Classic likelihood computation: full pair-hmm all haplotypes vs all reads.
 */
class PairHMMLikelihoodCalculationEngine {
private:
    const static int MAX_STR_UNIT_LENGTH = 8;
    const static int MAX_REPEAT_LENGTH   = 20;
    inline static int MIN_ADJUSTED_QSCORE = 10;


    double log10globalReadMismappingRate;
    std::shared_ptr<PairHMM> pairHMM;

    PCRErrorModel pcrErrorModel;

    char baseQualityScoreThreshold;

    /**
     * The expected rate of random sequencing errors for a read originating from its true haplotype.
     *
     * For example, if this is 0.01, then we'd expect 1 error per 100 bp.
     */
    static double EXPECTED_ERROR_RATE_PER_BASE;

    std::shared_ptr<char[]> pcrIndelErrorModelCache;

    void initializePCRErrorModel();

    static char getErrorModelAdjustedQual(int repeatLength, double rateFactor);

    /**
     * Initialize our pairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     *
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *
     * @param haplotypes a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    void initializePairHMM(std::vector<std::shared_ptr<Haplotype>> & haplotypes, std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> &perSampleReadList);

    void computeReadLikelihoods(SampleMatrix<SAMRecord, Haplotype> * likelihoods);

    /**
     * Pre-processing of the reads to be evaluated at the current location from the current sample.
     * We apply the PCR Error Model, and cap the minimum base, insertion, and deletion qualities of each read.
     * Modified copies of reads are packed into a new list, while original reads are retained for downstream use
     *
     * @param reads The original list of unmodified reads
     * @return processedReads. A new list of reads, in the same order, whose qualities have been altered by PCR error model and minimal quality thresholding
     */
    shared_ptr<vector<shared_ptr<SAMRecord>>> modifyReadQualities(vector<shared_ptr<SAMRecord>>& reads);

    /**
     *
     * @param length        the length of readBases
     * @param readBases
     * @param readInsQuals
     * @param readDelQuals
     */
    void applyPCRErrorModel(int length, uint8_t* readBases, uint8_t* readInsQuals, uint8_t* readDelQuals);

    static double log10MinTrueLikelihood(shared_ptr<SAMRecord> read, double maximumErrorPerBase);

    static int findTandemRepeatUnits(uint8_t* readBases, int length, int offset);

    static void capMinimumReadQualities(SAMRecord& read, int readQualsLength, uint8_t* readQuals, uint8_t* readInsQuals, uint8_t* readDelQuals, char baseQualityScoreThreshold);

    static uint8_t setToFixedValueIfTooLow(uint8_t currentVal, uint8_t minQual, uint8_t fixedQual);

    /**
    * Creates a new GATKRead with the source read's header, read group and mate
    * information, but with the following fields set to user-supplied values:
    *  - Read Bases
    *  - Base Qualities
    *  - Base Insertion Qualities
    *  - Base Deletion Qualities
    *
    *  Cigar string is empty (not-null)
    *
    * Use this method if you want to create a new GATKRead based on
    * another GATKRead, but with modified bases and qualities
    *
    * @param read a read to copy the header from
    * @param readBases an array containing the new bases you wish use in place of the originals
    * @param baseQualities an array containing the new base qualities you wish use in place of the originals
    * @param baseInsertionQualities an array containing the new base insertion qaulities
    * @param baseDeletionQualities an array containing the new base deletion qualities
    * @return a read with modified bases and qualities, safe for the GATK
    */
    static shared_ptr<SAMRecord> createQualityModifiedRead(SAMRecord& read, int length, std::shared_ptr<uint8_t[]> readBases, std::shared_ptr<uint8_t[]> baseQualities, std::shared_ptr<uint8_t[]> baseInsertionQualities, std::shared_ptr<uint8_t[]> baseDeletionQualities);

    static phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>* buildGapContinuationPenalties(vector<shared_ptr<SAMRecord>>& reads, char gapPenalty);

public:
    static double INITIAL_QSCORE;

    static char constantGCP;


    /**
     * Create a new PairHMMLikelihoodCalculationEngine using provided parameters and hmm to do its calculations
     *
     * @param constantGCP the gap continuation penalty to use with the PairHMM
     * @param hmmType the type of the HMM to use
     * @param log10globalReadMismappingRate the global mismapping probability, in log10(prob) units.  A value of
     *                                      -3 means that the chance that a read doesn't actually belong at this
     *                                      location in the genome is 1 in 1000.  The effect of this parameter is
     *                                      to cap the maximum likelihood difference between the reference haplotype
     *                                      and the best alternative haplotype by -3 log units.  So if the best
     *                                      haplotype is at -10 and this parameter has a value of -3 then even if the
     *                                      reference haplotype gets a score of -100 from the pairhmm it will be
     *                                      assigned a likelihood of -13.
     * @param pcrErrorModel model to correct for PCR indel artifacts
     * @param baseQualityScoreThreshold Base qualities below this threshold will be reduced to the minimum usable base
     */
    PairHMMLikelihoodCalculationEngine(char constantGCP, PairHMMNativeArgumentCollection& args, double log10globalReadMismappingRate, PCRErrorModel pcrErrorModel, char baseQualityScoreThreshold);

    ~PairHMMLikelihoodCalculationEngine();

    bool hasRateFactor() { return pcrErrorModel != 0; }

    double getRateFactor() { return (double)pcrErrorModel; }

    /**
     * Calculates the likelihood of reads across many samples evaluated against haplotypes resulting from the
     * active region assembly process.
     *
     * @param assemblyResultSet the input assembly results.
     * @param samples the list of targeted samples.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @throws IllegalArgumentException if any parameter is {@code null}.
     *
     * @return never {@code null}, and with at least one entry for input sample (keys in {@code perSampleReadList}.
     *    The value maps can be potentially empty though.
     */
    AlleleLikelihoods<SAMRecord, Haplotype>* computeReadLikelihoods(AssemblyResultSet & assemblyResultSet, std::vector<std::string>& samples, std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>& perSampleReadList);


};


#endif //MUTECT2CPP_MASTER_PAIRHMMLIKELIHOODCALCULATIONENGINE_H
