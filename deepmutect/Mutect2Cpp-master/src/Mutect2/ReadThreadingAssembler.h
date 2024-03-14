//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
#define MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H

#include "AssemblyResult.h"
#include "haplotype/Haplotype.h"
#include "SimpleInterval.h"
#include "AssemblyResultSet.h"
#include "tools/haplotypecaller/ReadErrorCorrector.h"
#include "path/ChainPruner.h"
#include <fstream>

class ReadThreadingAssembler {
private:
	bool recoverDanglingBranches = true;
	int pruneFactor;
	int numPruningSamples;
	int minDanglingBranchLength = 0;
	int numBestHaplotypesPerGraph;
	bool recoverAllDanglingBranches = false;
	bool removePathsNotConnectedToRef = true;
	bool justReturnRawGraph = false;
	bool debugGraphTransformations = false;
	std::vector<int> kmerSizes;
	bool dontIncreaseKmerSizesForCycles;
	int MAX_KMER_ITERATIONS_TO_ATTEMPT = 6;
	bool allowNonUniqueKmersInRef;
	ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> *chainPruner;
	static const uint8_t DEFAULT_MIN_BASE_QUALITY_TO_USE = 10;
	static const int KMER_SIZE_ITERATION_INCREASE = 10;
	std::string refSequenceName = "ref";

	std::shared_ptr<AssemblyResult> getAssemblyResult(std::shared_ptr<Haplotype> &refHaplotype,
	                                                  const std::shared_ptr<ReadThreadingGraph> &rtgraph) const;

	static std::shared_ptr<AssemblyResult> cleanupSeqGraph(const std::shared_ptr<SeqGraph> &seqGraph);

	/*std::vector<std::shared_ptr<Haplotype>>
	findBestPaths(const std::vector<std::shared_ptr<SeqGraph>> &graphs, std::shared_ptr<Haplotype> &refHaplotype,
	              const std::shared_ptr<SimpleInterval> &refLoc,
	              const std::shared_ptr<SimpleInterval> &activeRegionWindow,
	              const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> &assemblyResultByGraph,
	              std::shared_ptr<AssemblyResultSet> &assemblyResultSet) const;*/

	void findBestPaths(const std::vector<std::shared_ptr<SeqGraph>> &graphs, std::shared_ptr<Haplotype> &refHaplotype,
	                   const std::shared_ptr<SimpleInterval> &refLoc,
	                   const std::shared_ptr<SimpleInterval> &activeRegionWindow,
	                   const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> &assemblyResultByGraph,
	                   std::shared_ptr<AssemblyResultSet> &assemblyResultSet) const;

	std::shared_ptr<AssemblyResult>
	createGraph(const std::vector<std::shared_ptr<SAMRecord>> &reads, std::shared_ptr<Haplotype> &refHaplotype,
	            int kmerSize, bool allowLowComplexityGraphs);

	static void addResult(std::vector<std::shared_ptr<AssemblyResult>> &results,
	                      const std::shared_ptr<AssemblyResult> &maybeNullResult);

	static int getMinKmerSize(std::shared_ptr<Haplotype> &refHaplotype, std::vector<int> candidateKmerSizes);

public:
	ReadThreadingAssembler(int pruneFactor, int numPruningSamples, int numBestHaplotypesPerGraph,
	                       bool dontIncreaseKmerSizesForCycles, bool allowNonUniqueKmersInRef,
	                       std::vector<int> kmerSizes);

	ReadThreadingAssembler(int maxAllowedPathsForReadThreadingAssembler, std::vector<int> kmerSizes,
	                       bool dontIncreaseKmerSizesForCycles, bool allowNonUniqueKmersInRef, int numPruningSamples,
	                       int pruneFactor, bool useAdaptivePruning, double initialErrorRateForPruning,
	                       double pruningLogOddsThreshold, int maxUnprunedVariants);

	~ReadThreadingAssembler();

	std::shared_ptr<AssemblyResultSet>
	runLocalAssembly(const std::shared_ptr<AssemblyRegion> &assemblyRegion, std::shared_ptr<Haplotype> &refHaplotype,
	                 const std::shared_ptr<uint8_t[]> &fullReferenceWithPadding, int refLength,
	                 const std::shared_ptr<SimpleInterval> &refLoc, ReadErrorCorrector *readErrorCorrector);

	std::vector<std::shared_ptr<AssemblyResult>>
	assemble(std::vector<std::shared_ptr<SAMRecord>> &reads, std::shared_ptr<Haplotype> &refHaplotype);

	void setMinDanglingBranchLength(int minDanglingBranchLength);

	void setRecoverDanglingBranches(bool recoverAllDanglingBranches);

	void setRecoverAllDanglingBranches(bool recoverAllDanglingBranches);

protected:
	uint8_t minBaseQualityToUseInAssembly = DEFAULT_MIN_BASE_QUALITY_TO_USE;
};


#endif //MUTECT2CPP_MASTER_READTHREADINGASSEMBLER_H
