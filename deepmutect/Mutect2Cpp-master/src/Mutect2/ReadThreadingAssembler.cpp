//
// Created by 梦想家xixi on 2021/11/15.
//

#include "ReadThreadingAssembler.h"
#include <memory>
#include <utility>
#include "AssemblyResultSet.h"
#include "graph/KBestHaplotypeFinder.h"
#include "graph/ReadThreadingGraph.h"
#include "read/CigarUtils.h"
#include "AdaptiveChainPruner.h"
#include "boost/dynamic_bitset.hpp"
#include "graph/utils/GraphObjectPool.h"


std::shared_ptr<AssemblyResult>
ReadThreadingAssembler::getAssemblyResult(std::shared_ptr<Haplotype> &refHaplotype,
                                          const std::shared_ptr<ReadThreadingGraph> &rtgraph) const {
	if (recoverDanglingBranches) {
		if (pruneFactor < 0)
			throw std::invalid_argument("pruneFactor must be non-negative");
		if (minDanglingBranchLength < 0)
			throw std::invalid_argument("minDanglingBranchLength must be non-negative");
		if (!rtgraph->ifAlreadyBuilt())
			throw std::invalid_argument("recoverDanglingTails requires the graph be already built");
		rtgraph->recoverDanglingTails(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
		//rtgraph->printGraphSize("recoverDanglingTails");
		rtgraph->recoverDanglingHeads(pruneFactor, minDanglingBranchLength, recoverAllDanglingBranches);
		//rtgraph->printGraphSize("recoverDanglingHeads");
	}
	if (removePathsNotConnectedToRef) {
		rtgraph->removePathsNotConnectedToRef(rtgraph->getVertexSet().size());
		//rtgraph->printGraphSize("removePathsNotConnectedToRef");
	}

	std::shared_ptr<SeqGraph> initialSeqGraph = rtgraph->toSequenceGraph();
	//initialSeqGraph->printGraphSize("initialSeqGraph");
	if (justReturnRawGraph) {
		return std::make_shared<AssemblyResult>(ASSEMBLED_SOME_VARIATION, initialSeqGraph, nullptr);
	}
	initialSeqGraph->cleanNonRefPaths();
	//initialSeqGraph->printGraphSize("cleanNonRefPaths");

	std::shared_ptr<AssemblyResult> cleaned = cleanupSeqGraph(initialSeqGraph);
	return std::make_shared<AssemblyResult>(cleaned->getStatus(), cleaned->getGraph(), rtgraph);
}

std::shared_ptr<AssemblyResult> ReadThreadingAssembler::cleanupSeqGraph(const std::shared_ptr<SeqGraph> &seqGraph) {
	seqGraph->zipLinearChains();
	//seqGraph->printGraphSize("zipLinearChains");
	seqGraph->removeSingletonOrphanVertices();
	//seqGraph->printGraphSize("removeSingletonOrphanVertices");
	seqGraph->removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
	//seqGraph->printGraphSize("removeVerticesNotConnectedToRefRegardlessOfEdgeDirection");
	seqGraph->simplifyGraph();
	if (seqGraph->getReferenceSinkVertex() == nullptr || seqGraph->getReferenceSourceVertex() == nullptr) {
		return std::make_shared<AssemblyResult>(JUST_ASSEMBLED_REFERENCE, seqGraph, nullptr);
	}
	seqGraph->removePathsNotConnectedToRef(seqGraph->getVertexSet().size());
	//seqGraph->printGraphSize("removePathsNotConnectedToRef");
	seqGraph->simplifyGraph();
	if (seqGraph->getVertexSet().size() == 1) {
		std::shared_ptr<SeqVertex> complete = *(seqGraph->getVertexSet().begin());
		std::shared_ptr<SeqVertex> dummy(new SeqVertex(nullptr, 0));
		seqGraph->addVertex(dummy);
		seqGraph->addEdge(complete, dummy, GraphObjectPool::createSeqEdge(true, 0));
		//seqGraph->printGraphSize("dummy");
	}
	return std::make_shared<AssemblyResult>(ASSEMBLED_SOME_VARIATION, seqGraph, nullptr);
}

std::shared_ptr<AssemblyResultSet>
ReadThreadingAssembler::runLocalAssembly(const std::shared_ptr<AssemblyRegion> &assemblyRegion,
                                         std::shared_ptr<Haplotype> &refHaplotype,
                                         const std::shared_ptr<uint8_t[]> &fullReferenceWithPadding, int refLength,
                                         const std::shared_ptr<SimpleInterval> &refLoc,
                                         ReadErrorCorrector *readErrorCorrector) {
	Mutect2Utils::validateArg(assemblyRegion.get(), "Assembly engine cannot be used with a null AssemblyRegion.");
	Mutect2Utils::validateArg(refHaplotype.get(), "Active region must have an extended location.");
	Mutect2Utils::validateArg(fullReferenceWithPadding.get(), "fullReferenceWithPadding");
	Mutect2Utils::validateArg(refLoc.get(), "refLoc");
	Mutect2Utils::validateArg(refLength == refLoc->size(), "Reference bases and reference loc must be the same size.");

	std::vector<std::shared_ptr<SAMRecord>> correctedReads;
	if (readErrorCorrector != nullptr) {
		//TODO::readErrorCorrector
		readErrorCorrector->addReadsToKmers(assemblyRegion->getReads());
	}
	correctedReads = assemblyRegion->getReads();
	std::vector<std::shared_ptr<SeqGraph>> nonRefGraphs;
	std::shared_ptr<AssemblyResultSet> resultSet = std::make_shared<AssemblyResultSet>();
	resultSet->setRegionForGenotyping(assemblyRegion);
	resultSet->setFullReferenceWithPadding(fullReferenceWithPadding, refLength);
	resultSet->setPaddedReferenceLoc(refLoc);

	const std::shared_ptr<SimpleInterval> activeRegionExtendedLocation = assemblyRegion->getExtendedSpan();
	refHaplotype->setGenomeLocation(activeRegionExtendedLocation);
	resultSet->add(refHaplotype);

	std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> assemblyResultByGraph;
	for (const auto &result: assemble(correctedReads, refHaplotype)) {
		if (result->getStatus() == ASSEMBLED_SOME_VARIATION) {
			assemblyResultByGraph.insert(std::make_pair(result->getGraph(), result));
			nonRefGraphs.emplace_back(result->getGraph());
		}
	}
	/*std::cout << "nonRefGraphs\t" << nonRefGraphs.size() << std::endl;
	std::cout << "refHaplotype\t" << refHaplotype->getLength() << " " << refHaplotype->getBaseString() << std::endl;
	std::cout << "activeRegionExtendedLocation\t" << activeRegionExtendedLocation->getContig() << " "
	          << activeRegionExtendedLocation->getStart() + 1 << " " << activeRegionExtendedLocation->getEnd() + 1
	          << std::endl;*/

	findBestPaths(nonRefGraphs, refHaplotype, refLoc, activeRegionExtendedLocation, assemblyResultByGraph, resultSet);
	//resultSet->printSortedHaplotypes();
	return resultSet;
}

int
ReadThreadingAssembler::getMinKmerSize(std::shared_ptr<Haplotype> &refHaplotype, std::vector<int> candidateKmerSizes) {
	//std::string s = refHaplotype->getBaseString();
	std::shared_ptr<uint8_t[]> s = refHaplotype->getBases();
	bool justACGT = true;
	int len = refHaplotype->getLength();
	int i, k = 0;

	//when candidateKmerSizes[k] <= 30, use Bit Operation (long long)
	uint8_t charToU8[len];
	for (i = 0; i < len; ++i) {
		uint8_t ch = s[i];
		if (ch == 'A') charToU8[i] = 0;
		else if (ch == 'C') charToU8[i] = 1;
		else if (ch == 'G') charToU8[i] = 2;
		else if (ch == 'T') charToU8[i] = 3;
		else {
			justACGT = false;
			break;
		}
	}

	if (justACGT) {
		while (candidateKmerSizes[k] <= 30) {
            phmap::flat_hash_set<long long> valueSet;
			valueSet.reserve(len - candidateKmerSizes[k] + 1);
			long long val = 0L, mask = (1L << (candidateKmerSizes[k] * 2)) - 1;
			for (i = 0; i < candidateKmerSizes[k]; ++i) val = (val << 2) | charToU8[i];
			valueSet.insert(val);
			for (i = candidateKmerSizes[k]; i < len; ++i) {
				val = ((val << 2) & mask) | charToU8[i];
				if (!valueSet.insert(val).second) {
					k++;
					break;
				}
			}
			if (i == len) return candidateKmerSizes[k];
		}

		//when candidateKmerSizes[k] > 30, use dynamic bitset
		uint8_t charToBitH[len], charToBitL[len];
		for (i = 0; i < len; ++i) {
			uint8_t ch = s[i];
			if (ch == 'A') charToBitH[i] = 0, charToBitL[i] = 0;
			else if (ch == 'C') charToBitH[i] = 0, charToBitL[i] = 1;
			else if (ch == 'G') charToBitH[i] = 1, charToBitL[i] = 0;
			else charToBitH[i] = 1, charToBitL[i] = 1;
		}

		int last_i = 0, j;
		while (k < candidateKmerSizes.size() - 1) {
			int bitSetSize = 2 * candidateKmerSizes[k];
			boost::dynamic_bitset<> s1(bitSetSize), s2(bitSetSize);
			for (j = last_i; j < last_i + candidateKmerSizes[k] - 1; ++j) {
				s1 <<= 2;
				s1[1] = charToBitH[j], s1[0] = charToBitL[j];
			}
			for (i = last_i + candidateKmerSizes[k] - 1; i < len; i++) {
				s1 <<= 2;
				s1[1] = charToBitH[i], s1[0] = charToBitL[i];
				//std::cout << "s1: " << s1 << std::endl;
				s2 = s1;
				for (j = i + 1; j < len; j++) { //s2[j] loop
					s2 <<= 2;
					s2[1] = charToBitH[j], s2[0] = charToBitL[j];
					//std::cout << "s2: " << s2 << std::endl;
					if (s1 == s2) { //match
						last_i = i - candidateKmerSizes[k] + 1;
						k++;
						if (last_i + candidateKmerSizes[k] >= len) return candidateKmerSizes[k];
						break;
					}
				}
				if (j != len) break; //match, solve next k
			}
			if (i == len) return candidateKmerSizes[k]; //not match, return
		}
		return candidateKmerSizes[k];
	}

	// not justACGT, use 4 bits to solve
	uint8_t charToBit1[len], charToBit2[len], charToBit3[len], charToBit4[len]; //1 is highest, 2 is lowest
	for (i = 0; i < len; ++i) {
		uint8_t ch = s[i];
		if (ch == 'A') charToBit1[i] = 0, charToBit2[i] = 0, charToBit3[i] = 0, charToBit4[i] = 0;
		else if (ch == 'C') charToBit1[i] = 0, charToBit2[i] = 0, charToBit3[i] = 0, charToBit4[i] = 1;
		else if (ch == 'G') charToBit1[i] = 0, charToBit2[i] = 0, charToBit3[i] = 1, charToBit4[i] = 0;
		else if (ch == 'T') charToBit1[i] = 0, charToBit2[i] = 0, charToBit3[i] = 1, charToBit4[i] = 1;
		else if (ch == 'R') charToBit1[i] = 0, charToBit2[i] = 1, charToBit3[i] = 0, charToBit4[i] = 0;
		else if (ch == 'Y') charToBit1[i] = 0, charToBit2[i] = 1, charToBit3[i] = 0, charToBit4[i] = 1;
		else if (ch == 'M') charToBit1[i] = 0, charToBit2[i] = 1, charToBit3[i] = 1, charToBit4[i] = 0;
		else if (ch == 'K') charToBit1[i] = 0, charToBit2[i] = 1, charToBit3[i] = 1, charToBit4[i] = 1;
		else if (ch == 'S') charToBit1[i] = 1, charToBit2[i] = 0, charToBit3[i] = 0, charToBit4[i] = 0;
		else if (ch == 'W') charToBit1[i] = 1, charToBit2[i] = 0, charToBit3[i] = 0, charToBit4[i] = 1;
		else if (ch == 'H') charToBit1[i] = 1, charToBit2[i] = 0, charToBit3[i] = 1, charToBit4[i] = 0;
		else if (ch == 'B') charToBit1[i] = 1, charToBit2[i] = 0, charToBit3[i] = 1, charToBit4[i] = 1;
		else if (ch == 'V') charToBit1[i] = 1, charToBit2[i] = 1, charToBit3[i] = 0, charToBit4[i] = 0;
		else if (ch == 'D') charToBit1[i] = 1, charToBit2[i] = 1, charToBit3[i] = 0, charToBit4[i] = 1;
		else if (ch == 'N') charToBit1[i] = 1, charToBit2[i] = 1, charToBit3[i] = 1, charToBit4[i] = 0;
		else charToBit1[i] = 1, charToBit2[i] = 1, charToBit3[i] = 1, charToBit4[i] = 1;
	}

	int last_i = 0, j;
	while (k < candidateKmerSizes.size() - 1) {
		int bitSetSize = 4 * candidateKmerSizes[k];
		boost::dynamic_bitset<> s1(bitSetSize), s2(bitSetSize);
		for (j = last_i; j < last_i + candidateKmerSizes[k] - 1; ++j) {
			s1 <<= 4;
			s1[3] = charToBit1[j], s1[2] = charToBit2[j], s1[1] = charToBit3[j], s1[0] = charToBit4[j];
		}
		for (i = last_i + candidateKmerSizes[k] - 1; i < len; i++) {
			s1 <<= 4;
			s1[3] = charToBit1[i], s1[2] = charToBit2[i], s1[1] = charToBit3[i], s1[0] = charToBit4[i];
			//std::cout << "s1: " << s1 << std::endl;
			s2 = s1;
			for (j = i + 1; j < len; j++) { //s2[j] loop
				s2 <<= 4;
				s2[3] = charToBit1[j], s2[2] = charToBit2[j], s2[1] = charToBit3[j], s2[0] = charToBit4[j];
				//std::cout << "s2: " << s2 << std::endl;
				if (s1 == s2) { //match
					last_i = i - candidateKmerSizes[k] + 1;
					k++;
					if (last_i + candidateKmerSizes[k] >= len) return candidateKmerSizes[k];
					break;
				}
			}
			if (j != len) break; //match, solve next k
		}
		if (i == len) return candidateKmerSizes[k]; //not match, return
	}
	return candidateKmerSizes[k];
}

std::vector<std::shared_ptr<AssemblyResult>>
ReadThreadingAssembler::assemble(std::vector<std::shared_ptr<SAMRecord>> &reads,
                                 std::shared_ptr<Haplotype> &refHaplotype) {
	std::vector<std::shared_ptr<AssemblyResult>> results;
	if (kmerSizes.empty()) return results;

	std::sort(kmerSizes.begin(), kmerSizes.end());
	std::vector<int> tmpKmerSizes(kmerSizes);
	for (int numIterations = 0; numIterations < MAX_KMER_ITERATIONS_TO_ATTEMPT; ++numIterations) {
		tmpKmerSizes.push_back(*(tmpKmerSizes.end() - 1) + KMER_SIZE_ITERATION_INCREASE);
	}
	int minKmerSize = getMinKmerSize(refHaplotype, tmpKmerSizes);
	//std::cout << "minKmerSize = " << minKmerSize << std::endl;

	for (int kmerSize: kmerSizes) {
		if (kmerSize < minKmerSize && !allowNonUniqueKmersInRef) continue;
		addResult(results, createGraph(reads, refHaplotype, kmerSize, dontIncreaseKmerSizesForCycles));
	}

	if (results.empty() && !dontIncreaseKmerSizesForCycles) {
		int numIterations = 1, kmerSize = *(kmerSizes.end() - 1) + KMER_SIZE_ITERATION_INCREASE;
		while (numIterations < MAX_KMER_ITERATIONS_TO_ATTEMPT) {
			if (kmerSize >= minKmerSize || allowNonUniqueKmersInRef) {
				addResult(results, createGraph(reads, refHaplotype, kmerSize, false));
				if (!results.empty()) break;
			}
			numIterations++, kmerSize += KMER_SIZE_ITERATION_INCREASE;
		}
		if (numIterations == MAX_KMER_ITERATIONS_TO_ATTEMPT && results.empty())
			addResult(results, createGraph(reads, refHaplotype, kmerSize, true));
	}
	/*std::cout << "----------results Size " << results.size() << "----------\n";
	for (const auto &result: results) {
		result->getThreadingGraph()->printGraphSize("");
		//result->getThreadingGraph()->outputDotFile("threadingGraph.dot");
		result->getGraph()->printGraphSize("");
		//result->getGraph()->outputDotFile("seqGraph.dot");
		std::cout << "----------\n";
	}
	std::cout << "----------------------------------\n";*/
	return results;
}

std::shared_ptr<AssemblyResult>
ReadThreadingAssembler::createGraph(const std::vector<std::shared_ptr<SAMRecord>> &reads,
                                    std::shared_ptr<Haplotype> &refHaplotype, int kmerSize,
                                    bool allowLowComplexityGraphs) {
	if (refHaplotype->getLength() < kmerSize)
		return std::make_shared<AssemblyResult>(FAILED, nullptr, nullptr);

	std::shared_ptr<ReadThreadingGraph> rtgraph = std::make_shared<ReadThreadingGraph>(kmerSize,
	                                                                                   debugGraphTransformations,
	                                                                                   minBaseQualityToUseInAssembly,
	                                                                                   numPruningSamples);
	rtgraph->setThreadingStartOnlyAtExistingVertex(!recoverDanglingBranches);
	rtgraph->addSequence(refSequenceName, refHaplotype->getBases(), refHaplotype->getLength(), true);
	rtgraph->reserveSpace((int) (refHaplotype->getLength() * 1.1));

	for (std::shared_ptr<SAMRecord> read: reads)
		rtgraph->addRead(read);

#ifdef SORT_MODE
	rtgraph->sortPendingBySequence();
#endif

	//rtgraph->printPendingInfo();
	rtgraph->buildGraphIfNecessary();
	//rtgraph->printGraphSize("");
	//rtgraph->outputDotFile("./graph1.dot");

	chainPruner->pruneLowWeightChains(rtgraph);
	//rtgraph->printGraphSize("");
	//rtgraph->outputDotFile("./graph2.dot");

	if (rtgraph->hasCycles()) {
		//std::cout << std::to_string(kmerSize) + " failed because hasCycles" + '\n';
		return nullptr;
	}

	if (!allowLowComplexityGraphs && rtgraph->isLowComplexity()) {
		//std::cout << std::to_string(kmerSize) + " failed because isLowComplexity" + '\n';
		return nullptr;
	}
	//std::cout << std::to_string(kmerSize) + " OK\n";
	return getAssemblyResult(refHaplotype, rtgraph);
}

void ReadThreadingAssembler::addResult(std::vector<std::shared_ptr<AssemblyResult>> &results,
                                       const std::shared_ptr<AssemblyResult> &maybeNullResult) {
	if (maybeNullResult != nullptr) {
		results.emplace_back(maybeNullResult);
	}
}

void ReadThreadingAssembler::findBestPaths(const std::vector<std::shared_ptr<SeqGraph>> &graphs,
                                           std::shared_ptr<Haplotype> &refHaplotype,
                                           const std::shared_ptr<SimpleInterval> &refLoc,
                                           const std::shared_ptr<SimpleInterval> &activeRegionWindow,
                                           const std::map<std::shared_ptr<SeqGraph>, std::shared_ptr<AssemblyResult>> &assemblyResultByGraph,
                                           std::shared_ptr<AssemblyResultSet> &assemblyResultSet) const {
    phmap::flat_hash_set<std::shared_ptr<Haplotype>, hash_Haplotype, equal_Haplotype> returnHaplotypes;
	int activeRegionStart = refHaplotype->getAlignmentStartHapwrtRef();
	//int failedCigars = 0;

	for (auto &graph: graphs) {
		std::shared_ptr<SeqVertex> source = graph->getReferenceSourceVertex();
		std::shared_ptr<SeqVertex> sink = graph->getReferenceSinkVertex();
		Mutect2Utils::validateArg(source != nullptr && sink != nullptr, "Both source and sink cannot be null");

		for (auto &kBestHaplotype: KBestHaplotypeFinder(graph, source, sink).findBestHaplotypes(
				numBestHaplotypesPerGraph)) {
			std::shared_ptr<Haplotype> h = kBestHaplotype->getHaplotype();
			if (returnHaplotypes.find(h) == returnHaplotypes.end()) {
				if (kBestHaplotype->getIsReference()) {
					refHaplotype->setScore(kBestHaplotype->getScore());
				}
				std::shared_ptr<Cigar> cigar = CigarUtils::calculateCigar(refHaplotype->getBases(),
				                                                          refHaplotype->getLength(), h->getBases(),
				                                                          h->getLength());
				if (cigar == nullptr)
					continue;   // couldn't produce a meaningful alignment of haplotype to reference, fail quietly

				/*std::cout << refHaplotype->getLength() << " " <<std::string((char *)refHaplotype->getBases().get()) << std::endl << h->getLength()  << " " <<std::string((char *)h->getBases().get()) << std::endl;
				for (const auto &ce: cigar->getCigarElements()){
					std::cout << ce.getLength() << CigarOperatorUtils::enumToCharacter(ce.getOperator());
				}
				std::cout << std::endl;
				for (int i = 1; i < cigar->getCigarElements().size(); ++i){
					if (cigar->getCigarElement(i).getLength() == cigar->getCigarElement(i-1).getLength() && cigar->getCigarElement(i).getLength()>10){
						if((cigar->getCigarElement(i).getOperator()==D && cigar->getCigarElement(i-1).getOperator()==I)
						||(cigar->getCigarElement(i).getOperator()==I && cigar->getCigarElement(i-1).getOperator()==D)){
							std::shared_ptr<Cigar> cigar1 = CigarUtils::calculateCigar(refHaplotype->getBases(),refHaplotype->getLength(), h->getBases(),h->getLength());
							std::cout<<"break\n";
						}
					}
				}*/
				h->setCigar(cigar);
				h->setAlignmentStartHapwrtRef(activeRegionStart);
				h->setGenomeLocation(activeRegionWindow);
				returnHaplotypes.insert(h);
				assemblyResultSet->add(h, assemblyResultByGraph.at(graph));
			}
		}
	}
}

ReadThreadingAssembler::ReadThreadingAssembler(int pruneFactor, int numPruningSamples, int numBestHaplotypesPerGraph,
                                               bool dontIncreaseKmerSizesForCycles,
                                               bool allowNonUniqueKmersInRef, std::vector<int> kmerSizes)
		: pruneFactor(pruneFactor), numPruningSamples(numPruningSamples),
		  numBestHaplotypesPerGraph(numBestHaplotypesPerGraph),
		  dontIncreaseKmerSizesForCycles(dontIncreaseKmerSizesForCycles),
		  allowNonUniqueKmersInRef(allowNonUniqueKmersInRef), kmerSizes(std::move(kmerSizes)) {
	chainPruner = new AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(0.001, 2.302585092994046, 100);
	setMinDanglingBranchLength(4);
}

void ReadThreadingAssembler::setMinDanglingBranchLength(int minDanglingBranchLength) {
	this->minDanglingBranchLength = minDanglingBranchLength;
}

ReadThreadingAssembler::ReadThreadingAssembler(int maxAllowedPathsForReadThreadingAssembler, std::vector<int> kmerSizes,
                                               bool dontIncreaseKmerSizesForCycles, bool allowNonUniqueKmersInRef,
                                               int numPruningSamples, int pruneFactor, bool useAdaptivePruning,
                                               double initialErrorRateForPruning, double pruningLogOddsThreshold,
                                               int maxUnprunedVariants) : kmerSizes(std::move(kmerSizes)),
                                                                          dontIncreaseKmerSizesForCycles(
		                                                                          dontIncreaseKmerSizesForCycles),
                                                                          allowNonUniqueKmersInRef(
		                                                                          allowNonUniqueKmersInRef),
                                                                          numPruningSamples(numPruningSamples),
                                                                          pruneFactor(pruneFactor),
                                                                          numBestHaplotypesPerGraph(
		                                                                          maxAllowedPathsForReadThreadingAssembler) {
	if (maxAllowedPathsForReadThreadingAssembler < 1)
		throw std::invalid_argument("numBestHaplotypesPerGraph should be >= 1");
	chainPruner = new AdaptiveChainPruner<MultiDeBruijnVertex, MultiSampleEdge>(initialErrorRateForPruning,
	                                                                            pruningLogOddsThreshold,
	                                                                            maxUnprunedVariants);
}

ReadThreadingAssembler::~ReadThreadingAssembler() {
	delete chainPruner;
}

void ReadThreadingAssembler::setRecoverDanglingBranches(bool recoverDanglingBranches) {
	this->recoverDanglingBranches = recoverDanglingBranches;
}

void ReadThreadingAssembler::setRecoverAllDanglingBranches(bool recoverAllDanglingBranches) {
	this->recoverAllDanglingBranches = recoverAllDanglingBranches;
	recoverDanglingBranches = true;
}
