//
// Created by 梦想家xixi on 2022/3/1.
//

#ifndef MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H
#define MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H

#include "PairHMM.h"
#include "haplotype/Haplotype.h"
#include "samtools/SAMRecord.h"
#include "HaplotypeDataHolder.h"
#include "AssemblyResultSet.h"
#include "haplotypecaller/PairHMMNativeArgumentCollection.h"
#include "trie/buildTreeUtils.h"

class VectorLoglessPairHMM : public PairHMM {
private:
	std::vector<HaplotypeDataHolder> mHaplotypeDataArray;
	phmap::flat_hash_map<std::shared_ptr<Haplotype>, int, hash_Haplotype, equal_Haplotype> haplotypeToHaplotypeListIdxMap;
	unsigned mHaplotypeDataArrayLength;
	trieNode *root = nullptr;

public:
	explicit VectorLoglessPairHMM(PairHMMNativeArgumentCollection &args);

	~VectorLoglessPairHMM() override;

	/**
	 * Create a VectorLoglessPairHMM
	 *
	 * @param implementation    which implementation to use (AVX or OMP)
	 * @param args              arguments to the native GKL implementation
	 */
	void initialize(const std::vector<std::shared_ptr<Haplotype>> &haplotypes,
	                const std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> &perSampleReadList,
	                int readMaxLength, int haplotypeMaxLength) override;

	// unique reads version updated by hlf
	void computeLog10Likelihoods(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
	                             vector<shared_ptr<SAMRecord>> &processedReads,
	                             phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) override;

	// all reads but using trie version by lxy
	void computeLog10Likelihoods_trie(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
	                                      vector<shared_ptr<SAMRecord>> &processedReads,
	                                      phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) override;

	// combine trie version and unique version, edited by hlf
	void computeLog10Likelihoods_trie_unique(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
	                                             vector<shared_ptr<SAMRecord>> &processedReads,
	                                             phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) override;
};


#endif //MUTECT2CPP_MASTER_VECTORLOGLESSPAIRHMM_H
