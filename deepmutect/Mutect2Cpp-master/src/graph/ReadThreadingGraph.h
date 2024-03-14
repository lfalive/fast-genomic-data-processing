//
// Created by 梦想家xixi on 2021/10/18.
//

#ifndef MUTECT2CPP_MASTER_READTHREADINGGRAPH_H
#define MUTECT2CPP_MASTER_READTHREADINGGRAPH_H

#include "samtools/SAMRecord.h"
#include "MultiDeBruijnVertex.h"
#include "MultiSampleEdge.h"
#include <string>
#include <utility>
#include <vector>
#include <deque>
#include <list>
#include "Kmer.h"
#include "BaseGraph/DirectedSpecifics.h"
#include "DanglingChainMergeHelper.h"
#include "SeqGraph.h"

typedef struct {
	std::string name;
	std::shared_ptr<uint8_t[]> sequence;
	int start;
	int stop;
	int count;
	bool isRef;
} SequenceForKmers;

enum TraversalDirection {
	downwards,
	upwards
};

class ReadThreadingGraph : public DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge> {
protected:
	int kmerSize;

	DanglingChainMergeHelper *
	generateCigarAgainstDownwardsReferencePath(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
	                                           int minDanglingBranchLength, bool recoverAll);

	DanglingChainMergeHelper *
	generateCigarAgainstUpwardsReferencePath(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
	                                         int minDanglingBranchLength, bool recoverAll);

	std::shared_ptr<uint8_t[]>
	getBasesForPath(const std::vector<std::shared_ptr<MultiDeBruijnVertex>> &path, int &length, bool expandSource);

private:
	int numPruningSamples;

	std::shared_ptr<Kmer> refSource;

    phmap::flat_hash_set<std::shared_ptr<Kmer>, hash_kmer, equal_kmer> nonUniqueKmers;

    phmap::flat_hash_map<std::shared_ptr<Kmer>, std::shared_ptr<MultiDeBruijnVertex>, hash_kmer, equal_kmer> uniqueKmers;

	const uint8_t minBaseQualityToUseInAssembly;

	bool debugGraphTransformations{};

	bool startThreadingOnlyAtExistingVertex{};

	bool increaseCountsThroughBranches = false;

	static const bool INCREASE_COUNTS_BACKWARDS = true;

	static const int MAX_CIGAR_COMPLEXITY = 3;

	int maxMismatchesInDanglingHead = -1;

private:
	inline static std::string ANONYMOUS_SAMPLE = "XXX_UNNAMED_XXX";
	inline static std::string pendingKeys[3] = {ANONYMOUS_SAMPLE, "normal", "tumor"};

	/**
	 * Determines whether a base can safely be used for assembly.
	 * Currently disallows Ns and/or those with low quality
	 *
	 * @param base  the base under consideration
	 * @param qual  the quality of that base
	 * @return true if the base can be used for assembly, false otherwise
	 */
	bool baseIsUsableForAssembly(uint8_t base, uint8_t qual) const;

	bool alreadyBuilt{};

	std::map<std::string, std::vector<SequenceForKmers>> pending;

	void
	addSequence(std::string seqName, std::string &sampleName, const std::shared_ptr<uint8_t[]> &sequence, int start,
	            int stop, int count, bool isRef);

	/**
	 * Get the collection of all sequences for kmers across all samples in no particular order
	 * @return non-null Collection
	 */
	//std::list<SequenceForKmers> getAllPendingSequences();

	/**
	 * Compute the smallest kmer size >= minKmerSize and <= maxKmerSize that has no non-unique kmers
	 * among all sequences added to the current graph.  Will always return a result for maxKmerSize if
	 * all smaller kmers had non-unique kmers.
	 *
	 * @param minKmerSize the minimum kmer size to consider when constructing the graph
	 * @param maxKmerSize the maximum kmer size to consider
	 * @return a non-null NonUniqueResult
	 */

	/**
	* Create a new vertex for kmer.  Add it to the uniqueKmers map if appropriate.
	*
	* kmer must not have a entry in unique kmers, or an error will be thrown
	*
	* @param kmer the kmer we want to create a vertex for
	* @return the non-null created vertex
	*/
	std::shared_ptr<MultiDeBruijnVertex> createVertex(const std::shared_ptr<Kmer> &kmer);

	/**
   * Workhorse routine of the assembler.  Given a sequence whose last vertex is anchored in the graph, extend
   * the graph one bp according to the bases in sequence.
   *
   * @param prevVertex a non-null vertex where sequence was last anchored in the graph
   * @param sequence the sequence we're threading through the graph
   * @param kmerStart the start of the current kmer in graph we'd like to add
   * @param count the number of observations of this kmer in graph (can be > 1 for GGA)
   * @param isRef is this the reference sequence?
   * @return a non-null vertex connecting prevVertex to in the graph based on sequence
   */
	std::shared_ptr<MultiDeBruijnVertex>
	extendChainByOne(const std::shared_ptr<MultiDeBruijnVertex> &prevVertex, const std::shared_ptr<uint8_t[]> &sequence,
	                 int kmerStart, int count, bool isRef);

	void threadSequence(SequenceForKmers &sequenceForKmers);

	/**
	* Find vertex and its position in seqForKmers where we should start assembling seqForKmers
	*
	* @param seqForKmers the sequence we want to thread into the graph
	* @return the position of the starting vertex in seqForKmer, or -1 if it cannot find one
	*/
	int findStart(const SequenceForKmers &seqForKmers);

	std::shared_ptr<MultiDeBruijnVertex> getUniqueKmerVertex(const std::shared_ptr<Kmer> &kmer, bool allowRefSource);

	std::shared_ptr<MultiDeBruijnVertex> getOrCreateKmerVertex(const std::shared_ptr<uint8_t[]> &sequence, int start);

	void increaseCountsInMatchedKmers(int incr, const std::shared_ptr<MultiDeBruijnVertex> &vertex,
	                                  const std::shared_ptr<uint8_t[]> &originalKmer, int offset);

	/**
	* Attempt to attach vertex with out-degree == 0 to the graph
	*
	* @param vertex the vertex to recover
	* @param pruneFactor  the prune factor to use in ignoring chain pieces
	* @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
	* @param aligner
	* @return 1 if we successfully recovered the vertex and 0 otherwise
	*/
	/*int recoverDanglingTail(const std::shared_ptr<MultiDeBruijnVertex> &v, int pruneFactor, int minDanglingBranchLength,
	                        bool recoverAll);*/


	std::vector<std::shared_ptr<MultiDeBruijnVertex>>
	findPathUpwardsToLowestCommonAncestor(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
	                                      bool giveUpAtBranch);

	bool hasIncidentRefEdge(const std::shared_ptr<MultiDeBruijnVertex> &v);

	std::shared_ptr<MultiSampleEdge> getHeaviestIncomingEdge(const std::shared_ptr<MultiDeBruijnVertex> &v);

	std::shared_ptr<MultiSampleEdge> getHeaviestOutgoingEdge(const std::shared_ptr<MultiDeBruijnVertex> &v);

	std::vector<std::shared_ptr<MultiDeBruijnVertex>>
	getReferencePath(std::shared_ptr<MultiDeBruijnVertex> start, TraversalDirection direction,
	                 const std::shared_ptr<MultiSampleEdge> &blacklistedEdge);

	/**
	 * Attempt to attach vertex with in-degree == 0, or a vertex on its path, to the graph
	 *
	 * @param vertex the vertex to recover
	 * @param pruneFactor  the prune factor to use in ignoring chain pieces
	 * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
	 * @param recoverAll recover even branches with forks
	 * @param aligner
	 * @return 1 if we successfully recovered a vertex and 0 otherwise
	 */
	/*int recoverDanglingHead(const std::shared_ptr<MultiDeBruijnVertex> &v, int pruneFactor, int minDanglingBranchLength,
	                        bool recoverAll);*/

	std::vector<std::shared_ptr<MultiDeBruijnVertex>>
	findPathDownwardsToHighestCommonDescendantOfReference(std::shared_ptr<MultiDeBruijnVertex> vertex, int pruneFactor,
	                                                      bool giveUpAtBranch);

	int bestPrefixMatch(const std::shared_ptr<uint8_t[]> &path1, const std::shared_ptr<uint8_t[]> &path2, int maxIndex);

	int getMaxMismatches(int lengthOfDanglingBranch) const;

	bool extendDanglingPathAgainstReference(DanglingChainMergeHelper *danglingHeadMergeResult, int numNodesToExtend);

	void resetToInitialState();

public:
	ReadThreadingGraph(uint8_t minBaseQualityToUseInAssembly, int kmerSize, bool alreadyBuilt,
	                   std::shared_ptr<Kmer> ref, int numPruningSamples)
			: minBaseQualityToUseInAssembly(minBaseQualityToUseInAssembly), kmerSize(kmerSize),
			  alreadyBuilt(false), refSource(std::move(ref)), numPruningSamples(numPruningSamples), nonUniqueKmers() {}

	ReadThreadingGraph(int kmerSize, bool debugGraphTransformations, uint8_t minBaseQualityToUseInAssembly,
	                   int numPruningSamples);

	/**
	 * Build the read threaded assembly graph if it hasn't already been constructed from the sequences that have
	 * been added to the graph.
	 */
	void buildGraphIfNecessary();

	bool ifAlreadyBuilt() const;

	/**
	 * Changes the threading start location policy.
	 *
	 * @param value  {@code true} if threading will start only at existing vertices in the graph, {@code false} if
	 *  it can start at any unique kmer.
	 */
	void setThreadingStartOnlyAtExistingVertex(bool value);

	/**
	 * Add a read to the sequence graph.  Finds maximal consecutive runs of bases with sufficient quality
	 * and applies {@see addSequence} to these subreads if they are longer than the kmer size.
	 *
	 * @param read a non-null read
	 */
	void addRead(std::shared_ptr<SAMRecord> &read);

	bool removeVertex(const std::shared_ptr<MultiDeBruijnVertex> &V) override;

	//void setPending();


	void removeSingletonOrphanVertices() override;

	/**
	 * Try to recover dangling tails
	 *
	 * @param pruneFactor  the prune factor to use in ignoring chain pieces
	 * @param minDanglingBranchLength the minimum length of a dangling branch for us to try to merge it
	 * @param recoverAll recover even branches with forks
	 * @param aligner
	 */
	void recoverDanglingTails(int pruneFactor, int minDanglingBranchLength, bool recoverAll);

	static bool cigarIsOkayToMerge(std::shared_ptr<Cigar> &cigar, bool requireFirstElementM, bool requireLastElementM);

	void recoverDanglingHeads(int pruneFactor, int minDanglingBranchLength, bool recoverAll);

	/**
	 * Actually merge the dangling tail if possible
	 *
	 * @param danglingTailMergeResult   the result from generating a Cigar for the dangling tail against the reference
	 * @return 1 if merge was successful, 0 otherwise
	 */
	int mergeDanglingTail(DanglingChainMergeHelper *danglingTailMergeResult);

	static int
	longestSuffixMatch(const std::shared_ptr<uint8_t[]> &seq, int seqLength, const std::shared_ptr<uint8_t[]> &kmer,
	                   int kmerLength, int seqStart);

	int mergeDanglingHead(DanglingChainMergeHelper *danglingTailMergeResult);

	std::shared_ptr<MultiSampleEdge>
	createEdge(const std::shared_ptr<MultiDeBruijnVertex> &, const std::shared_ptr<MultiDeBruijnVertex> &) override;

	std::shared_ptr<SeqGraph> toSequenceGraph();

	void
	addSequence(std::string seqName, const std::shared_ptr<uint8_t[]> &sequence, int length, int count, bool isRef);

	void addSequence(std::string seqName, const std::shared_ptr<uint8_t[]> &sequence, int length, bool isRef);

	bool isLowComplexity();

	bool hasCycles();

	void determineNonUniques();

	void reserveSpace(int size) override;

	void printPendingInfo();

	void printGraphSize(const std::string &info);

	void sortPendingBySequence();
};


#endif //MUTECT2CPP_MASTER_READTHREADINGGRAPH_H
