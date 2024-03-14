//
// Class for performing the pair HMM for global alignment. Figure 4.1 in Durbin 1998 book.
// Created by lhh on 4/22/22.
//

#ifndef MUTECT2CPP_MASTER_PAIRHMM_H
#define MUTECT2CPP_MASTER_PAIRHMM_H

#include "Haplotype.h"
#include "samtools/SAMRecord.h"
#include "utils/genotyper/AlleleLikelihoods.h"

class PairHMM {
public:
    PairHMM() = default;
    virtual ~PairHMM();

    const static char BASE_QUALITY_SCORE_THRESHOLD = 18;    // Base quals less than this value are squashed down to min possible qual

    virtual void initialize(const std::vector<std::shared_ptr<Haplotype>> & haplotypes,
                    const std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> & perSampleReadList,
                    int readMaxLength, int haplotypeMaxLength);

    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     *
     * Note: Do not worry about padding, just provide the true max length of the read and haplotype. The HMM will take care of the padding.
     *
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     * @throws IllegalArgumentException if haplotypeMaxLength is less than or equal to zero
     */
    void initialize(int readMaxLength, int haplotypeMaxLength);

    /**
     *  Given a list of reads and haplotypes, for every read compute the total probability of said read arising from
     *  each haplotype given base substitution, insertion, and deletion probabilities.
     *
     * @param processedReads reads to analyze instead of the ones present in the destination read-likelihoods.
     * @param logLikelihoods where to store the log likelihoods where position [a][r] is reserved for the log likelihood of {@code reads[r]}
     *             conditional to {@code alleles[a]}.
     * @param gcp penalty for gap continuations base array map for processed reads.
     *
     */
    virtual void computeLog10Likelihoods(SampleMatrix<SAMRecord, Haplotype>* logLikelihoods,
                                 vector<shared_ptr<SAMRecord>>& processedReads,
                                         phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>* gcp) = 0;

    virtual void computeLog10Likelihoods_trie(SampleMatrix<SAMRecord, Haplotype>* logLikelihoods,
                                         vector<shared_ptr<SAMRecord>>& processedReads,
                                         phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>* gcp) = 0;

	virtual void computeLog10Likelihoods_trie_unique(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
	                                                     vector<shared_ptr<SAMRecord>> &processedReads,
	                                                     phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) = 0;
    bool is_use_trietree_optimize;

protected:
    bool constantsAreInitialized = false;
    int hapStartIndex;

    int maxHaplotypeLength, maxReadLength;
    int paddedMaxReadLength, paddedMaxHaplotypeLength;
    int paddedReadLength, paddedHaplotypeLength;
    bool initialized = false;

    vector<double> mLogLikelihoodArray_1D;
	std::vector<std::vector<double>> mLogLikelihoodArray_2D;
};


#endif //MUTECT2CPP_MASTER_PAIRHMM_H
