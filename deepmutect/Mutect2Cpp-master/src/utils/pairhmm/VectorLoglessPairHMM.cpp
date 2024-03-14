//
// Created by 梦想家xixi on 2022/3/1.
//

#include <iomanip>
#include "VectorLoglessPairHMM.h"
#include "intel/pairhmm/IntelPairHmm.h"
#include "ReadUtils.h"
#include "haplotypecaller/ReadForPairHMM.h"
#include "parallel_hashmap/phmap.h"
#include "utils/pairhmm/PairHMMConcurrentControl.h"
#include "intel/common/avx.h"

VectorLoglessPairHMM::VectorLoglessPairHMM(PairHMMNativeArgumentCollection &args) : mHaplotypeDataArrayLength(0) {
	initNative(args.useDoublePrecision, args.pairHmmNativeThreads);
}

void VectorLoglessPairHMM::initialize(const std::vector<std::shared_ptr<Haplotype>> &haplotypes,
                                      const std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> &perSampleReadList,
                                      const int readMaxLength, const int haplotypeMaxLength) {

	mHaplotypeDataArrayLength = haplotypes.size();
	mHaplotypeDataArray.clear();
	mHaplotypeDataArray.reserve(mHaplotypeDataArrayLength);
	haplotypeToHaplotypeListIdxMap.clear();
	haplotypeToHaplotypeListIdxMap.reserve(mHaplotypeDataArrayLength);
	int idx = 0;
	phmap::flat_hash_set<int> length_count;
	for (const std::shared_ptr<Haplotype> &currHaplotype: haplotypes) {
		int len = currHaplotype->getBasesLength();
		length_count.insert(len);
		mHaplotypeDataArray.emplace_back(currHaplotype->getBases().get(), len);
		haplotypeToHaplotypeListIdxMap.emplace(currHaplotype, idx++);
	}
	is_use_trietree_optimize = haplotypes.size() / length_count.size() > 3;

	if (is_use_trietree_optimize) {
		buildTreeUtils::deleteTree(root);
		root = buildTreeUtils::buildTreeWithHaplotype_same_height(haplotypes, true);
	}
}

void VectorLoglessPairHMM::computeLog10Likelihoods(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
                                                   vector<shared_ptr<SAMRecord>> &processedReads,
                                                   phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) {
	if (processedReads.empty())
		return;

	int numReads = processedReads.size();

	// An array that stores unique testcases and where all testcases are mapped to this array
	// Therefore, the number of testcases to be computed is reduced, but the length of mloglikelihoodarray is unchanged
	std::vector<testcase> uniqueTestcases;
	/* For Debugging*/
//	std::vector<testcase> old_uniqueTestcases;
	std::vector<int> mapAlltoUnique;
	uniqueTestcases.reserve(numReads * mHaplotypeDataArrayLength);
	mapAlltoUnique.reserve(numReads * mHaplotypeDataArrayLength);

	// Where do all the testcases related to a read start
	phmap::flat_hash_map<std::shared_ptr<ReadForPairHMM>, int, ReadForPairHMMHash, ReadForPairHMMEqual> uniqueReadForPairHMM;
	uniqueReadForPairHMM.reserve(numReads);

	// Generate unique testcases
	int _rslen;
	char *gapConts;
	uint8_t *reads, *readQuals;
	std::vector<shared_ptr<uint8_t[]>> insGops, delGops;
	insGops.reserve(numReads);
	delGops.reserve(numReads);
	for (int r = 0; r < numReads; ++r) {
		_rslen = processedReads[r]->getLength();
		insGops.emplace_back(ReadUtils::getBaseInsertionQualities(processedReads[r], _rslen));
		delGops.emplace_back(ReadUtils::getBaseDeletionQualities(processedReads[r], _rslen));
		gapConts = (*gcp)[processedReads[r].get()].get();
		readQuals = processedReads[r]->getBaseQualitiesNoCopy().get();
		reads = processedReads[r]->getBasesNoCopy().get();

		std::shared_ptr<ReadForPairHMM> readOfTestcase = std::make_shared<ReadForPairHMM>(_rslen, readQuals,
		                                                                                  insGops[r].get(),
		                                                                                  delGops[r].get(), gapConts,
		                                                                                  reads);
		auto readIt = uniqueReadForPairHMM.find(readOfTestcase);
		if (BOOST_LIKELY(readIt == uniqueReadForPairHMM.end())) {
			// Push testcases into uniqueTestcases and mark the index where the first testcase appears
			uniqueReadForPairHMM.emplace(readOfTestcase, uniqueTestcases.size());
			readOfTestcase->initializeFloatVector();    // initialize probaility arrays and distm
			for (int h = 0; h < mHaplotypeDataArrayLength; ++h) {
				mapAlltoUnique.emplace_back(uniqueTestcases.size());
				uniqueTestcases.emplace_back(mHaplotypeDataArray[h].length, mHaplotypeDataArray[h].haplotypeBases,
				                             readOfTestcase);
				/* For Debugging*/
//				old_uniqueTestcases.emplace_back(haplotypeLengths[h], haplotypes[h], readOfTestcase);
			}
		} else {
			// No element needs to be pushed into uniqueTestcases
			int mapStart = readIt->second;
			for (int h = 0; h < mHaplotypeDataArrayLength; ++h) {
				mapAlltoUnique.emplace_back(mapStart++);
				/* For Debugging*/
//				old_uniqueTestcases.emplace_back(haplotypeLengths[h], haplotypes[h], readOfTestcase);
			}
		}
	}

//	PairHMMConcurrentControl::unique_reads += uniqueReadForPairHMM.size();
//	PairHMMConcurrentControl::all_reads += numReads;
//	PairHMMConcurrentControl::unique_cases += uniqueTestcases.size();
//	PairHMMConcurrentControl::all_cases += numReads * mHaplotypeDataArrayLength;
//	if (uniqueTestcases.size() != numReads * mHaplotypeDataArrayLength) {
//		std::cout << "==========================\n";
//		std::cout << "read:\t" << uniqueReadForPairHMM.size() << " / " << numReads << std::endl;
//		std::cout << "case:\t" << uniqueTestcases.size() << " / " << numReads * mHaplotypeDataArrayLength << std::endl;
//	}

	// Compute
	std::vector<double> uniqueLogLikelihoodArray(uniqueTestcases.size());
	for (int i = 0; i < uniqueTestcases.size(); ++i)
		computeLikelihoodsNative_concurrent_i(uniqueTestcases, uniqueLogLikelihoodArray, i);
//	computeLikelihoodsNative_concurrent(uniqueTestcases, uniqueLogLikelihoodArray);

	// Mapping results
	mLogLikelihoodArray_1D.resize(mapAlltoUnique.size());
	for (int i = 0; i < mapAlltoUnique.size(); ++i)
		mLogLikelihoodArray_1D[i] = uniqueLogLikelihoodArray[mapAlltoUnique[i]];

	/* For Debugging*/
//	std::vector<double> old_mLogLikelihoodArray(old_uniqueTestcases.size());
//	computeLikelihoodsNative_concurrent(old_uniqueTestcases, old_mLogLikelihoodArray);
//	assert(old_mLogLikelihoodArray == mLogLikelihoodArray);

	std::cout.setf(ios::fixed);

	int readIdx = 0;
	for (int r = 0; r < numReads; ++r) {
		int hapIdx = 0;
		for (auto &haplotype: logLikelihoods->alleles()) {
			//Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
			//get idx of current haplotype in the list and use this idx to get the right likelihoodValue
			int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.at(haplotype);
			//std::cout << setprecision(5) << mLogLikelihoodArray_1D[readIdx + idxInsideHaplotypeList] << " ";
			logLikelihoods->set(hapIdx, r, mLogLikelihoodArray_1D[readIdx + idxInsideHaplotypeList]);
			hapIdx++;
		}
		readIdx += mHaplotypeDataArrayLength;
	}
	//std::cout << std::endl;
}

void VectorLoglessPairHMM::computeLog10Likelihoods_trie(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
                                                        vector<shared_ptr<SAMRecord>> &processedReads,
                                                        phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) {
	if (processedReads.empty())
		return;

	int numReads = processedReads.size();

	std::vector<trie_testcase> uniqueTestcases;
	uniqueTestcases.reserve(numReads * mHaplotypeDataArrayLength);

	// Generate unique testcases
	int _rslen;
	char *gapConts;
	uint8_t *reads, *readQuals;
	std::vector<shared_ptr<uint8_t[]>> insGops, delGops;
	insGops.reserve(numReads);
	delGops.reserve(numReads);
	for (int r = 0; r < numReads; ++r) {
		_rslen = processedReads[r]->getLength();
		insGops.emplace_back(ReadUtils::getBaseInsertionQualities(processedReads[r], _rslen));
		delGops.emplace_back(ReadUtils::getBaseDeletionQualities(processedReads[r], _rslen));
		gapConts = (*gcp)[processedReads[r].get()].get();
		readQuals = processedReads[r]->getBaseQualitiesNoCopy().get();
		reads = processedReads[r]->getBasesNoCopy().get();

		std::shared_ptr<ReadForPairHMM> readOfTestcase = std::make_shared<ReadForPairHMM>(_rslen, readQuals,
		                                                                                  insGops[r].get(),
		                                                                                  delGops[r].get(), gapConts,
		                                                                                  reads);
		readOfTestcase->initializeFloatVector();
		uniqueTestcases.emplace_back(mHaplotypeDataArray, readOfTestcase, root);
	}

	// Compute
	mLogLikelihoodArray_2D.clear();
	mLogLikelihoodArray_2D.resize(uniqueTestcases.size());
	for (auto &array: mLogLikelihoodArray_2D) {
		array.reserve(mHaplotypeDataArrayLength);
	}
	computeLikelihoodsNative_concurrent_trie(uniqueTestcases, mLogLikelihoodArray_2D);

	//std::cout.setf(ios::fixed);

	for (int r = 0; r < numReads; ++r) {
		int hapIdx = 0;
		for (auto &haplotype: logLikelihoods->alleles()) {
			//Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
			//get idx of current haplotype in the list and use this idx to get the right likelihoodValue
			int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.at(haplotype);
			//std::cout << setprecision(5) << mLogLikelihoodArray_2D[r][idxInsideHaplotypeList] << " ";
			logLikelihoods->set(hapIdx++, r, mLogLikelihoodArray_2D[r][idxInsideHaplotypeList]);
		}
	}
	//std::cout << std::endl;
}

void VectorLoglessPairHMM::computeLog10Likelihoods_trie_unique(SampleMatrix<SAMRecord, Haplotype> *logLikelihoods,
                                                               vector<shared_ptr<SAMRecord>> &processedReads,
                                                               phmap::flat_hash_map<SAMRecord *, shared_ptr<char[]>> *gcp) {
	if (processedReads.empty())
		return;

	int numReads = processedReads.size();

	std::vector<trie_testcase> uniqueTestcases;
	uniqueTestcases.reserve(numReads);
	std::vector<int> mapAlltoUnique;
	mapAlltoUnique.reserve(numReads);

	// map the results of unique trie_testcases to all trie_testcases
	phmap::flat_hash_map<std::shared_ptr<ReadForPairHMM>, int, ReadForPairHMMHash, ReadForPairHMMEqual> uniqueReadForPairHMM;
	uniqueReadForPairHMM.reserve(numReads);

	// Generate unique testcases
	int _rslen;
	char *gapConts;
	uint8_t *reads, *readQuals;
	std::vector<shared_ptr<uint8_t[]>> insGops, delGops;
	insGops.reserve(numReads);
	delGops.reserve(numReads);
	for (int r = 0; r < numReads; ++r) {
		_rslen = processedReads[r]->getLength();
		insGops.emplace_back(ReadUtils::getBaseInsertionQualities(processedReads[r], _rslen));
		delGops.emplace_back(ReadUtils::getBaseDeletionQualities(processedReads[r], _rslen));
		gapConts = (*gcp)[processedReads[r].get()].get();
		readQuals = processedReads[r]->getBaseQualitiesNoCopy().get();
		reads = processedReads[r]->getBasesNoCopy().get();

		std::shared_ptr<ReadForPairHMM> readOfTestcase = std::make_shared<ReadForPairHMM>(_rslen, readQuals,
		                                                                                  insGops[r].get(),
		                                                                                  delGops[r].get(), gapConts,
		                                                                                  reads);
		auto readIt = uniqueReadForPairHMM.find(readOfTestcase);
		if (BOOST_LIKELY(readIt == uniqueReadForPairHMM.end())) {
			// Push testcases into uniqueTestcases and mark the index where the testcase appears
			uniqueReadForPairHMM.emplace(readOfTestcase, uniqueTestcases.size());
			mapAlltoUnique.emplace_back(uniqueTestcases.size());
			uniqueTestcases.emplace_back(mHaplotypeDataArray, readOfTestcase, root);
			readOfTestcase->initializeFloatVector();    // initialize probaility arrays and distm
		} else {
			// No element needs to be pushed into uniqueTestcases
			mapAlltoUnique.emplace_back(readIt->second);
		}
	}

	// Count the number of nodes to be calculated by trie and the number of nodes in the original version.
//	PairHMMConcurrentControl::trie_nodes += (buildTreeUtils::numberOfNodes(root) - 1) * numReads;
//	int cnt = 0;
//	for (int h = 0; h < mHaplotypeDataArrayLength; ++h) {
//		cnt += mHaplotypeDataArray[h].length / 32 + 1;
//	}
//	PairHMMConcurrentControl::original_nodes += cnt * numReads;

	// Compute
	mLogLikelihoodArray_2D.clear();
	mLogLikelihoodArray_2D.resize(uniqueTestcases.size());
	for (auto &array: mLogLikelihoodArray_2D) {
		array.reserve(mHaplotypeDataArrayLength);
	}
	computeLikelihoodsNative_concurrent_trie(uniqueTestcases, mLogLikelihoodArray_2D);

	//std::cout.setf(ios::fixed);

	for (int r = 0; r < numReads; ++r) {
		int target = mapAlltoUnique[r];
		int hapIdx = 0;
		for (auto &haplotype: logLikelihoods->alleles()) {
			//Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
			//get idx of current haplotype in the list and use this idx to get the right likelihoodValue
			int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.at(haplotype);
			//std::cout << setprecision(5) << mLogLikelihoodArray_2D[target][idxInsideHaplotypeList] << " ";
			logLikelihoods->set(hapIdx++, r, mLogLikelihoodArray_2D[target][idxInsideHaplotypeList]);
		}
	}
	//std::cout << std::endl;
}

VectorLoglessPairHMM::~VectorLoglessPairHMM() {
	buildTreeUtils::deleteTree(root);
}


