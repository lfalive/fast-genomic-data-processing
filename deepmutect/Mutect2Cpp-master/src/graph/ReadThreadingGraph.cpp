//
// Created by 梦想家xixi on 2021/10/18.
//

#include "ReadThreadingGraph.h"
#include <memory>
#include <string>
#include <utility>
#include "Mutect2Utils.h"
#include "smithwaterman/SmithWatermanAlignment.h"
#include "smithwaterman/SWNativeAlignerWrapper.h"
#include "read/AlignmentUtils.h"
#include "SeqVertex.h"
#include "boost/dynamic_bitset.hpp"
#include "graph/utils/GraphObjectPool.h"

void ReadThreadingGraph::addRead(std::shared_ptr<SAMRecord> &read) {
	std::shared_ptr<uint8_t[]> sequence_ = read->getBasesNoCopy();
	uint8_t *sequence = sequence_.get();
	uint8_t *qualities = read->getBaseQualitiesNoCopy().get();

//    for(int i = 0; i < read->getLength(); i++) {
//        std::cout << (char)sequence[i];
//    }
//    std::cout << std::endl;
//
//    if(sequence[0] == 'A' && sequence[1] == 'G' && sequence[2] == 'G' && sequence[3] == 'A')
//        std::cout << "hello";

	int lastGood = -1, length = read->getLength();
	for (int end = 0; end <= length; end++) {
		if (end == length || !baseIsUsableForAssembly(sequence[end], qualities[end])) {
			int start = lastGood, len = end - start;
			if (start != -1 && len >= kmerSize) {
				std::string name = read->getName() + '_' + std::to_string(start) + '_' + std::to_string(end);
				std::string sampleName = read->getGroup() == 0 ? "normal" : "tumor";
				addSequence(name, sampleName, sequence_, start, end, 1, false);
			}
			lastGood = -1;
		} else if (lastGood == -1) {
			lastGood = end;
		}
	}
}

bool ReadThreadingGraph::baseIsUsableForAssembly(uint8_t base, uint8_t qual) const {
	return base != 'N' && qual >= minBaseQualityToUseInAssembly;
}

void ReadThreadingGraph::addSequence(std::string seqName, std::string &sampleName,
                                     const std::shared_ptr<uint8_t[]> &sequence, int start, int stop,
                                     int count, bool isRef) {
	if (alreadyBuilt)
		throw std::invalid_argument("Graph already built");
	auto iter = pending.find(sampleName);
	if (iter == pending.end()) {
		std::vector<SequenceForKmers> list;
		pending.insert(std::make_pair(sampleName, list));
		iter = pending.find(sampleName);
	}
	iter->second.push_back(SequenceForKmers{std::move(seqName), sequence, start, stop, count, isRef});
}

void ReadThreadingGraph::reserveSpace(int size) {
	DirectedSpecifics::reserveSpace(size);
	uniqueKmers.reserve(size);
}

//find kmers repeated in the same sequence
void ReadThreadingGraph::determineNonUniques() {
	nonUniqueKmers.clear();
	std::vector<std::pair<std::shared_ptr<uint8_t[]>, int>> nonUniqueKmersPair;

	// check justACGT
	bool justACGT = true;
	for (auto &iter: pending) {
		for (auto &withNonUnique: iter.second) {
			uint8_t *s = withNonUnique.sequence.get();
			int len = withNonUnique.stop - withNonUnique.start;
			for (int i = 0; i < len; ++i) {
				uint8_t ch = s[withNonUnique.start + i];
				if (BOOST_UNLIKELY(ch != 'G' && ch != 'T' && ch != 'A' && ch != 'C')) {
					justACGT = false;
					break;
				}
			}
			if (!justACGT) break;
		}
		if (!justACGT) break;
	}
	//std::cout<<"justACGT "+std::to_string(justACGT)+'\n';

	if (justACGT) {
		if (kmerSize <= 30) {  //when kmerSize_ <= 30, use Bit Operation (long long)
			phmap::flat_hash_set<long long> nonUniquesFromSeqSet;
			phmap::flat_hash_set<long long> seqAllKmers;

			for (auto &iter: pending) {
				for (auto &withNonUnique: iter.second) {
					int len = withNonUnique.stop - withNonUnique.start;
					seqAllKmers.clear();
					if (seqAllKmers.bucket_count() < len - kmerSize + 1) {
						seqAllKmers.reserve(len - kmerSize + 1);
					}

					//std::string str = reinterpret_cast<const char *>(withNonUnique.sequence.get());
					//sequence[start, stop)
					uint8_t *s = withNonUnique.sequence.get();
					uint8_t charToU8[len];
					for (int i = 0; i < len; ++i) {
						uint8_t ch = s[withNonUnique.start + i];
						if (ch == 'A') charToU8[i] = 0;
						else if (ch == 'C') charToU8[i] = 1;
						else if (ch == 'G') charToU8[i] = 2;
						else if (ch == 'T') charToU8[i] = 3;
						else throw std::invalid_argument("Found N in sequence when determining NonUniques!");
					}

					long long val = 0L, mask = (1L << (kmerSize * 2)) - 1;
					for (int i = 0; i < kmerSize - 1; ++i) val = (val << 2) | charToU8[i];
					for (int i = kmerSize - 1; i < len; i++) {
						val = ((val << 2) & mask) | charToU8[i];
						//std::cout << val << std::endl;
						if (!seqAllKmers.insert(val).second) {   // is nonUnique in the sequence
							if (nonUniquesFromSeqSet.insert(val).second) {   //and first appearance
								nonUniqueKmersPair.emplace_back(withNonUnique.sequence,
								                                withNonUnique.start + i - kmerSize + 1);
								//std::cout << withNonUnique.start << " " << withNonUnique.stop << " "
								//          << str.substr(withNonUnique.start, len) << std::endl;
								//std::cout << str.substr(withNonUnique.start + i - kmerSize_ + 1, kmerSize_) << std::endl;
							}
						}
					}
				}
			}
		} else {    //when kmerSize_ > 30, use dynamic bitset
			phmap::flat_hash_set<boost::dynamic_bitset<>> nonUniquesFromSeqSet;
			phmap::flat_hash_set<boost::dynamic_bitset<>> seqAllKmers;

			for (auto &iter: pending) {
				for (auto &withNonUnique: iter.second) {
					int len = withNonUnique.stop - withNonUnique.start;
					seqAllKmers.clear();
					if (seqAllKmers.bucket_count() < len - kmerSize + 1) {
						seqAllKmers.reserve(len - kmerSize + 1);
					}

					//std::string str = reinterpret_cast<const char *>(withNonUnique.sequence.get());
					uint8_t *s = withNonUnique.sequence.get();
					uint8_t charToBitH[len], charToBitL[len];
					for (int i = 0; i < len; ++i) {
						uint8_t ch = s[withNonUnique.start + i];
						if (ch == 'A') charToBitH[i] = 0, charToBitL[i] = 0;
						else if (ch == 'C') charToBitH[i] = 0, charToBitL[i] = 1;
						else if (ch == 'G') charToBitH[i] = 1, charToBitL[i] = 0;
						else if (ch == 'T') charToBitH[i] = 1, charToBitL[i] = 1;
						else throw std::invalid_argument("Found N in sequence when determining NonUniques!");
					}

					boost::dynamic_bitset<> bitSeq(2 * kmerSize);
					for (int i = 0; i < kmerSize - 1; ++i) {
						bitSeq <<= 2;
						bitSeq[1] = charToBitH[i], bitSeq[0] = charToBitL[i];
					}
					for (int i = kmerSize - 1; i < len; i++) {
						bitSeq <<= 2;
						bitSeq[1] = charToBitH[i], bitSeq[0] = charToBitL[i];
						//std::cout << bitSeq << std::endl;
						if (!seqAllKmers.insert(bitSeq).second) {   // is nonUnique in the sequence
							if (nonUniquesFromSeqSet.insert(bitSeq).second) {   //and first appearance
								nonUniqueKmersPair.emplace_back(withNonUnique.sequence,
								                                withNonUnique.start + i - kmerSize + 1);
								//std::cout << withNonUnique.start << " " << withNonUnique.stop << " "
								//          << str.substr(withNonUnique.start, len) << std::endl;
								//std::cout << str.substr(withNonUnique.start + i - kmerSize_ + 1, kmerSize_) << std::endl;
							}
						}
					}
				}
			}
		}
	} else {   // not justACGT, use 4 bits to solve
		phmap::flat_hash_set<boost::dynamic_bitset<>> nonUniquesFromSeqSet;
		phmap::flat_hash_set<boost::dynamic_bitset<>> seqAllKmers;

		for (auto &iter: pending) {
			for (auto &withNonUnique: iter.second) {
				int len = withNonUnique.stop - withNonUnique.start;
				seqAllKmers.clear();
				if (seqAllKmers.bucket_count() < len - kmerSize + 1) {
					seqAllKmers.reserve(len - kmerSize + 1);
				}

				//std::string str = reinterpret_cast<const char *>(withNonUnique.sequence.get());
				uint8_t *s = withNonUnique.sequence.get();
				uint8_t charToBit1[len], charToBit2[len], charToBit3[len], charToBit4[len]; //1 is highest, 2 is lowest
				for (int i = 0; i < len; ++i) {
					uint8_t ch = s[withNonUnique.start + i];
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

				boost::dynamic_bitset<> bitSeq(4 * kmerSize);
				for (int i = 0; i < kmerSize - 1; ++i) {
					bitSeq <<= 4;
					bitSeq[3] = charToBit1[i], bitSeq[2] = charToBit2[i], bitSeq[1] = charToBit3[i], bitSeq[0] = charToBit4[i];
				}
				for (int i = kmerSize - 1; i < len; i++) {
					bitSeq <<= 4;
					bitSeq[3] = charToBit1[i], bitSeq[2] = charToBit2[i], bitSeq[1] = charToBit3[i], bitSeq[0] = charToBit4[i];
					//std::cout << bitSeq << std::endl;
					if (!seqAllKmers.insert(bitSeq).second) {   // is nonUnique in the sequence
						if (nonUniquesFromSeqSet.insert(bitSeq).second) {   //and first appearance
							nonUniqueKmersPair.emplace_back(withNonUnique.sequence,
							                                withNonUnique.start + i - kmerSize + 1);
							//std::cout << withNonUnique.start << " " << withNonUnique.stop << " "
							//          << str.substr(withNonUnique.start, len) << std::endl;
							//std::cout << str.substr(withNonUnique.start + i - kmerSize_ + 1, kmerSize_) << std::endl;
						}
					}
				}
			}
		}
	}

	for (auto &p: nonUniqueKmersPair) {
		nonUniqueKmers.insert(std::make_shared<Kmer>(p.first, p.second, kmerSize));
	}
}

std::shared_ptr<MultiDeBruijnVertex> ReadThreadingGraph::createVertex(const std::shared_ptr<Kmer> &kmer) {
	std::shared_ptr<MultiDeBruijnVertex> newVertex = GraphObjectPool::createMultiVertex(kmer->getBases(),
	                                                                                   kmer->getLength(), false);

	unsigned prevSize = getVertexSet().size();
	addVertex(newVertex);
	if (getVertexSet().size() != prevSize + 1)
		throw std::invalid_argument("Adding vertex to graph didn't increase the graph size");

	if (nonUniqueKmers.find(kmer) == nonUniqueKmers.end())
		uniqueKmers.insert({kmer, newVertex});

	return newVertex;
}

std::shared_ptr<MultiDeBruijnVertex>
ReadThreadingGraph::extendChainByOne(const std::shared_ptr<MultiDeBruijnVertex> &prevVertex,
                                     const std::shared_ptr<uint8_t[]> &sequence, const int kmerStart, const int count,
                                     const bool isRef) {
	std::shared_ptr<MultiDeBruijnVertex> nextVertex;

	for (const auto &outgoingEdge: outgoingEdgesOf(prevVertex)) {
		nextVertex = getEdgeTarget(outgoingEdge);
		if (nextVertex->getSuffix() == sequence.get()[kmerStart + kmerSize - 1]) {
			outgoingEdge->incMultiplicity(count);
			return nextVertex;
		}
	}

	std::shared_ptr<Kmer> kmer = std::make_shared<Kmer>(sequence, kmerStart, kmerSize);
	if (*kmer == *refSource && !isRef) {
		nextVertex = createVertex(kmer);
	} else {
		nextVertex = getUniqueKmerVertex(kmer, false);
		if (nextVertex == nullptr) {
			nextVertex = createVertex(kmer);
		} else if (isRef)
			throw std::invalid_argument("Found a unique vertex to merge into the reference graph");
	}

	addEdge(prevVertex, nextVertex, GraphObjectPool::createMultiEdge(isRef, count, numPruningSamples));
	return nextVertex;
}

void ReadThreadingGraph::threadSequence(SequenceForKmers &sequenceForKmers) {
	int uniqueStartPos = findStart(sequenceForKmers);
	if (uniqueStartPos == -1)
		return;

	std::shared_ptr<MultiDeBruijnVertex> vertex = getOrCreateKmerVertex(sequenceForKmers.sequence,
	                                                                    uniqueStartPos);

	if (INCREASE_COUNTS_BACKWARDS)
		increaseCountsInMatchedKmers(sequenceForKmers.count, vertex, vertex->getSequence(), kmerSize - 2);

	if (sequenceForKmers.isRef) {
		if (refSource->getBases() != nullptr)
			throw std::invalid_argument("Found two refSources! prev:");
		// update refSource
		refSource = std::make_shared<Kmer>(sequenceForKmers.sequence, sequenceForKmers.start, kmerSize);
	}

	for (int i = uniqueStartPos + 1; i <= sequenceForKmers.stop - kmerSize; i++) {
		vertex = extendChainByOne(vertex, sequenceForKmers.sequence, i, sequenceForKmers.count,
		                          sequenceForKmers.isRef);
	}
}

int ReadThreadingGraph::findStart(const SequenceForKmers &seqForKmers) {
	if (seqForKmers.isRef)
		return 0;
	for (int i = seqForKmers.start; i < seqForKmers.stop - kmerSize; i++) {
		std::shared_ptr<Kmer> kmer1 = std::make_shared<Kmer>(seqForKmers.sequence, i, kmerSize);
		// when startThreadingOnlyAtExistingVertex, find kmer in uniqueKmers set
		// when !startThreadingOnlyAtExistingVertex, not find kmer in nonUniqueKmers set
		if ((startThreadingOnlyAtExistingVertex ? uniqueKmers.find(kmer1) != uniqueKmers.end() :
		     nonUniqueKmers.find(kmer1) == nonUniqueKmers.end()))
			return i;
	}
	return -1;
}

std::shared_ptr<MultiDeBruijnVertex>
ReadThreadingGraph::getUniqueKmerVertex(const std::shared_ptr<Kmer> &kmer, bool allowRefSource) {
	if (!allowRefSource && kmer == refSource)
		return nullptr;
	auto res = uniqueKmers.find(kmer);
	return res == uniqueKmers.end() ? nullptr : res->second;
}

std::shared_ptr<MultiDeBruijnVertex>
ReadThreadingGraph::getOrCreateKmerVertex(const std::shared_ptr<uint8_t[]> &sequence, int start) {
	std::shared_ptr<Kmer> kmer = std::make_shared<Kmer>(sequence, start, kmerSize);
	std::shared_ptr<MultiDeBruijnVertex> res = getUniqueKmerVertex(kmer, true);
	return res == nullptr ? createVertex(kmer) : res;
}

void
ReadThreadingGraph::increaseCountsInMatchedKmers(int incr, const std::shared_ptr<MultiDeBruijnVertex> &vertex,
                                                 const std::shared_ptr<uint8_t[]> &originalKmer, int offset) {
	//if (offset == -1)
	//	return;

//    if(vertex->getHashCode() == 884547439)
//        std::cout << "hello";

	std::queue<std::pair<std::shared_ptr<MultiDeBruijnVertex>, int>> q;
	q.emplace(vertex, offset);

	while (!q.empty()) {
		std::shared_ptr<MultiDeBruijnVertex> v = q.front().first;
		int o = q.front().second;
		q.pop();
		for (const auto &incomingEdge: incomingEdgesOf(v)) {
			std::shared_ptr<MultiDeBruijnVertex> prev = getEdgeSource(incomingEdge);
			if (prev->getSuffix() == originalKmer.get()[o] &&
			    (inDegreeOf(v) == 1 || increaseCountsThroughBranches)) {
				incomingEdge->incMultiplicity(incr);
				if (o - 1 >= 0)
					q.emplace(prev, o - 1);
			}
		}
	}
}

void ReadThreadingGraph::buildGraphIfNecessary() {
	if (alreadyBuilt)
		return;

	determineNonUniques();

#ifdef SORT_MODE
	for (const auto &key: pendingKeys) {
		if (pending.find(key) != pending.end()) {
			for (auto &viter: pending[key]) {
				threadSequence(viter);
				//outputDotFile("./CPP_" + std::to_string(cnt) + ".dot");
				//printGraphSize();
			}
			for (auto &eiter: edgeMap) {
				eiter.first->flushSingleSampleMultiplicity();
			}
		}
	}
#else
	for (auto &item: pending) {
		for (auto &viter: item.second) {
			threadSequence(viter);
		}
		for (auto &eiter: edgeMap) {
			eiter.first->flushSingleSampleMultiplicity();
		}
	}
#endif

	pending.clear();
	alreadyBuilt = true;

	for (auto &uniqueKmer: uniqueKmers) {
		uniqueKmer.second->additionalInfoAppendPlusSign();
	}
}

void ReadThreadingGraph::setThreadingStartOnlyAtExistingVertex(bool value) {
	startThreadingOnlyAtExistingVertex = value;
}

bool ReadThreadingGraph::removeVertex(const std::shared_ptr<MultiDeBruijnVertex> &V) {
	std::shared_ptr<uint8_t[]> sequence(new uint8_t[V->getLength()]);
	memcpy(sequence.get(), V->getSequence().get(), V->getLength());
	bool result = DirectedSpecifics::removeVertex(V);

	if (result) {
		const std::shared_ptr<Kmer> kmer = std::make_shared<Kmer>(sequence, 0, kmerSize);
		uniqueKmers.erase(kmer);
	}

	return result;
}

void ReadThreadingGraph::removeSingletonOrphanVertices() {
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> toRemove;
	phmap::flat_hash_set<std::shared_ptr<MultiDeBruijnVertex>> &allvertex = getVertexSet();
	typename phmap::flat_hash_set<std::shared_ptr<MultiDeBruijnVertex>>::iterator viter;
	for (viter = allvertex.begin(); viter != allvertex.end(); viter++) {
		if (inDegreeOf(*viter) == 0 && outDegreeOf(*viter) == 0) {
			toRemove.emplace_back(*viter);
		}
	}
	removeAllVertices(toRemove);
}

void ReadThreadingGraph::recoverDanglingTails(int pruneFactor, int minDanglingBranchLength, bool recoverAll) {
	DanglingChainMergeHelper *danglingTailMergeResult = nullptr;

#ifdef SORT_MODE
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> vertexSet = getSortedVertexList();
#else
	phmap::flat_hash_set<std::shared_ptr<MultiDeBruijnVertex>> vertexSet = getVertexSet();
#endif
	for (const std::shared_ptr<MultiDeBruijnVertex> &v: vertexSet) {
		if (outDegreeOf(v) == 0 && !isRefSink(v)) {
			danglingTailMergeResult = generateCigarAgainstDownwardsReferencePath(v, pruneFactor,
			                                                                     minDanglingBranchLength,
			                                                                     recoverAll);
			if (danglingTailMergeResult != nullptr) {
				if (cigarIsOkayToMerge(danglingTailMergeResult->cigar, false, true))
					mergeDanglingTail(danglingTailMergeResult);
				delete danglingTailMergeResult;
			}
		}
	}
}

DanglingChainMergeHelper *
ReadThreadingGraph::generateCigarAgainstDownwardsReferencePath(std::shared_ptr<MultiDeBruijnVertex> vertex,
                                                               int pruneFactor,
                                                               int minDanglingBranchLength, bool recoverAll) {
	int minTailPathLength = std::max(1, minDanglingBranchLength);
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> altPath = findPathUpwardsToLowestCommonAncestor(
			std::move(vertex), pruneFactor, !recoverAll);
	if (altPath.empty() || isRefSource(altPath.back()) || altPath.size() < minTailPathLength + 1)
		return nullptr;
	std::reverse(altPath.begin(), altPath.end());
	std::shared_ptr<MultiDeBruijnVertex> toBlacklistedEdge = altPath.size() > 1 ? altPath[1] : nullptr;
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> refPath = getReferencePath(altPath[0], downwards,
	                                                                             getHeaviestIncomingEdge(
			                                                                             toBlacklistedEdge));
	int refLength, altLength;
	std::shared_ptr<uint8_t[]> refBases = getBasesForPath(refPath, refLength, false);
	std::shared_ptr<uint8_t[]> altBases = getBasesForPath(altPath, altLength, false);
	SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
	SmithWatermanAlignment *alignment = wrapper.align(refBases, refLength, altBases, altLength,
	                                                  &SmithWatermanAligner::STANDARD_NGS, LEADING_INDEL);
	std::shared_ptr<Cigar> c = alignment->getCigar();
	delete alignment;
	return new DanglingChainMergeHelper(altPath, refPath, altBases, altLength, refBases, refLength,
	                                    AlignmentUtils::removeTrailingDeletions(c));
}

std::vector<std::shared_ptr<MultiDeBruijnVertex>>
ReadThreadingGraph::findPathUpwardsToLowestCommonAncestor(std::shared_ptr<MultiDeBruijnVertex> vertex,
                                                          int pruneFactor,
                                                          bool giveUpAtBranch) {
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> ret;
	std::shared_ptr<MultiDeBruijnVertex> v = std::move(vertex);
	if (giveUpAtBranch) {
		while (!(inDegreeOf(v) != 1 || outDegreeOf(v) >= 2)) {
			std::shared_ptr<MultiSampleEdge> edge = incomingEdgeOf(v);
			if (edge->getPruningMultiplicity() < pruneFactor) {
				ret.clear();
			} else {
				ret.push_back(v);
			}
			v = getEdgeSource(edge);
		}
		ret.push_back(v);

		return outDegreeOf(v) > 1 ? ret : std::vector<std::shared_ptr<MultiDeBruijnVertex>>();
	} else {
		while (!(hasIncidentRefEdge(v) || inDegreeOf(v) == 0)) {
			std::shared_ptr<MultiSampleEdge> edge = getHeaviestIncomingEdge(v);
			if (edge->getPruningMultiplicity() < pruneFactor) {
				ret.clear();
			} else {
				ret.push_back(v);
			}
			v = getEdgeSource(edge);
		}
		ret.push_back(v);

		return outDegreeOf(v) > 1 && hasIncidentRefEdge(v) ? ret : std::vector<std::shared_ptr<MultiDeBruijnVertex>>();
	}
}

bool ReadThreadingGraph::hasIncidentRefEdge(const std::shared_ptr<MultiDeBruijnVertex> &v) {
	return (std::any_of(incomingEdgesOf(v).begin(), incomingEdgesOf(v).end(),
	                    [](const std::shared_ptr<MultiSampleEdge> &edge) { return edge->getIsRef(); }));
	/*for (const std::shared_ptr<MultiSampleEdge> &edge: incomingEdgesOf(v)) {
		if (edge->getIsRef()) {
			return true;
		}
	}
	return false;*/
}

std::shared_ptr<MultiSampleEdge>
ReadThreadingGraph::getHeaviestIncomingEdge(const std::shared_ptr<MultiDeBruijnVertex> &v) {
	phmap::flat_hash_set<std::shared_ptr<MultiSampleEdge>> incomingEdges = incomingEdgesOf(v);
	std::shared_ptr<MultiSampleEdge> ret;
	ret = *incomingEdges.begin();
	for (const auto &incomingEdge: incomingEdges) {
		if (ret->getPruningMultiplicity() < incomingEdge->getPruningMultiplicity())
			ret = incomingEdge;
	}
	return ret;
}

std::vector<std::shared_ptr<MultiDeBruijnVertex>>
ReadThreadingGraph::getReferencePath(std::shared_ptr<MultiDeBruijnVertex> start, TraversalDirection direction,
                                     const std::shared_ptr<MultiSampleEdge> &blacklistedEdge) {
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> path;
	path.reserve(vertexMapDirected.size());
	std::shared_ptr<MultiDeBruijnVertex> v = std::move(start);
	while (v != nullptr) {
		path.emplace_back(v);
		v = (direction == downwards ?
		     getNextReferenceVertex(v, true, blacklistedEdge) : getPrevReferenceVertex(v));
	}
	return path;
}

std::shared_ptr<uint8_t[]>
ReadThreadingGraph::getBasesForPath(const std::vector<std::shared_ptr<MultiDeBruijnVertex>> &path, int &length,
                                    bool expandSource) {
	int tmpLength = 300;
	int start = 0;
	std::shared_ptr<uint8_t[]> tmp(new uint8_t[300]);
	for (const std::shared_ptr<MultiDeBruijnVertex> &v: path) {
		if (expandSource && isSource(v)) {
			std::shared_ptr<uint8_t[]> seq = v->getSequence();
			int seqLength = v->getLength();
			while (start + seqLength > tmpLength) {
				std::shared_ptr<uint8_t[]> newtmp(new uint8_t[2 * tmpLength]);
				memcpy(newtmp.get(), tmp.get(), tmpLength);
				tmpLength *= 2;
				tmp = newtmp;
			}
			for (int i = 0; i < seqLength; i++) {
				tmp.get()[start + i] = seq.get()[seqLength - 1 - i];
			}
			start += seqLength;
		} else {
			while (start + 1 > tmpLength) {
				int oldLength = tmpLength;
				tmpLength *= 2;
				std::shared_ptr<uint8_t[]> newtmp(new uint8_t[tmpLength]);
				memcpy(newtmp.get(), tmp.get(), oldLength);
				tmp = newtmp;
			}
			tmp.get()[start] = v->getSuffix();
			start++;
		}
	}
	std::shared_ptr<uint8_t[]> ret(new uint8_t[start + 1]);
	memcpy(ret.get(), tmp.get(), start);
	ret.get()[start] = '\0';
	length = start;
	return ret;
}

bool ReadThreadingGraph::cigarIsOkayToMerge(std::shared_ptr<Cigar> &cigar, bool requireFirstElementM,
                                            bool requireLastElementM) {
	std::vector<CigarElement> elements = cigar->getCigarElements();
	int numElements = elements.size();
	if (numElements == 0 || numElements > MAX_CIGAR_COMPLEXITY) {
		return false;
	}

	if (requireFirstElementM && elements.at(0).getOperator() != M) {
		return false;
	}

	if (requireLastElementM && elements.at(numElements - 1).getOperator() != M) {
		return false;
	}

	return true;
}

int ReadThreadingGraph::mergeDanglingTail(DanglingChainMergeHelper *danglingTailMergeResult) {
	std::vector<CigarElement> elements = danglingTailMergeResult->cigar->getCigarElements();
	CigarElement lastElement = elements.at(elements.size() - 1);
	Mutect2Utils::validateArg(lastElement.getOperator() == M, "The last Cigar element must be an M");

	int lastRefIndex = danglingTailMergeResult->cigar->getReferenceLength() - 1;
	int matchingSuffix = std::min(longestSuffixMatch(danglingTailMergeResult->referencePathString,
	                                                 danglingTailMergeResult->referencePathStringLength,
	                                                 danglingTailMergeResult->danglingPathString,
	                                                 danglingTailMergeResult->danglingPathStringLength,
	                                                 lastRefIndex), lastElement.getLength());

	if (matchingSuffix == 0) {
		return 0;
	}

	int altIndexToMerge = std::max(danglingTailMergeResult->cigar->getReadLength() - matchingSuffix - 1, 0);

	bool firstElementIsDeletion = elements.at(0).getOperator() == D;
	bool mustHandleLeadingDeletionCase =
			firstElementIsDeletion && (elements.at(0).getLength() + matchingSuffix == lastRefIndex + 1);
	int refIndexToMerge = lastRefIndex - matchingSuffix + 1 + (mustHandleLeadingDeletionCase ? 1 : 0);

	if (refIndexToMerge == 0) {
		return 0;
	}
	addEdge(danglingTailMergeResult->danglingPath.at(altIndexToMerge),
	        danglingTailMergeResult->referencePath.at(refIndexToMerge), GraphObjectPool::createMultiEdge(
					false, 1, numPruningSamples));
	return 1;
}

int ReadThreadingGraph::longestSuffixMatch(const std::shared_ptr<uint8_t[]> &seq, int seqLength,
                                           const std::shared_ptr<uint8_t[]> &kmer, int kmerLength,
                                           int seqStart) {
	for (int len = 1; len <= kmerLength; len++) {
		int seqI = seqStart - len + 1;
		int kmerI = kmerLength - len;
		if (seqI < 0 || seq.get()[seqI] != kmer.get()[kmerI]) {
			return len - 1;
		}
	}
	return kmerLength;
}

void ReadThreadingGraph::recoverDanglingHeads(int pruneFactor, int minDanglingBranchLength, bool recoverAll) {
	DanglingChainMergeHelper *danglingHeadMergeResult = nullptr;

#ifdef SORT_MODE
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> vertexSet = getSortedVertexList();
#else
	phmap::flat_hash_set<std::shared_ptr<MultiDeBruijnVertex>> vertexSet = getVertexSet();
#endif
	for (const std::shared_ptr<MultiDeBruijnVertex> &v: vertexSet) {
		if (inDegreeOf(v) == 0 && !isRefSource(v)) {
			danglingHeadMergeResult = generateCigarAgainstUpwardsReferencePath(v, pruneFactor,
			                                                                   minDanglingBranchLength,
			                                                                   recoverAll);
			if (danglingHeadMergeResult != nullptr) {
				if (cigarIsOkayToMerge(danglingHeadMergeResult->cigar, true, false))
					mergeDanglingHead(danglingHeadMergeResult);
				delete danglingHeadMergeResult;
			}
		}
	}
}

DanglingChainMergeHelper *
ReadThreadingGraph::generateCigarAgainstUpwardsReferencePath(std::shared_ptr<MultiDeBruijnVertex> vertex,
                                                             int pruneFactor,
                                                             int minDanglingBranchLength, bool recoverAll) {
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> altPath = findPathDownwardsToHighestCommonDescendantOfReference(
			std::move(vertex), pruneFactor, !recoverAll);
	if (altPath.empty() || isRefSink(altPath.back()) || altPath.size() < minDanglingBranchLength + 1) {
		return nullptr;
	}
	std::reverse(altPath.begin(), altPath.end());
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> refPath = getReferencePath(altPath[0], upwards, nullptr);
	int refLength, altLength;
	std::shared_ptr<uint8_t[]> refBases = getBasesForPath(refPath, refLength, true);
	std::shared_ptr<uint8_t[]> altBases = getBasesForPath(altPath, altLength, true);
	SWNativeAlignerWrapper wrapper = SWNativeAlignerWrapper();
	SmithWatermanAlignment *alignment = wrapper.align(refBases, refLength, altBases, altLength,
	                                                  &SmithWatermanAligner::STANDARD_NGS, LEADING_INDEL);
	std::shared_ptr<Cigar> c = alignment->getCigar();
	delete alignment;
	return new DanglingChainMergeHelper(altPath, refPath, altBases, altLength, refBases, refLength,
	                                    AlignmentUtils::removeTrailingDeletions(c));
}

std::vector<std::shared_ptr<MultiDeBruijnVertex>>
ReadThreadingGraph::findPathDownwardsToHighestCommonDescendantOfReference(std::shared_ptr<MultiDeBruijnVertex> vertex,
                                                                          int pruneFactor, bool giveUpAtBranch) {
	std::vector<std::shared_ptr<MultiDeBruijnVertex>> ret;
	std::shared_ptr<MultiDeBruijnVertex> v = std::move(vertex);
	if (giveUpAtBranch) {
		while (!(isReferenceNode(v) || outDegreeOf(v) != 1)) {
			std::shared_ptr<MultiSampleEdge> edge = outgoingEdgeOf(v);
			if (edge->getPruningMultiplicity() < pruneFactor) {
				ret.clear();
			} else {
				ret.push_back(v);
			}
			v = getEdgeTarget(edge);
		}
		ret.push_back(v);

		return isReferenceNode(v) ? ret : std::vector<std::shared_ptr<MultiDeBruijnVertex>>();
	} else {
		while (!(isReferenceNode(v) || outDegreeOf(v) == 0)) {
			std::shared_ptr<MultiSampleEdge> edge = getHeaviestOutgoingEdge(v);
			if (edge->getPruningMultiplicity() < pruneFactor) {
				ret.clear();
			} else {
				ret.push_back(v);
			}
			v = getEdgeTarget(edge);
		}
		ret.push_back(v);

		return isReferenceNode(v) ? ret : std::vector<std::shared_ptr<MultiDeBruijnVertex>>();
	}
}

std::shared_ptr<MultiSampleEdge>
ReadThreadingGraph::getHeaviestOutgoingEdge(const std::shared_ptr<MultiDeBruijnVertex> &v) {
	phmap::flat_hash_set<std::shared_ptr<MultiSampleEdge>> outgoing = outgoingEdgesOf(v);
	std::shared_ptr<MultiSampleEdge> ret;
	ret = *outgoing.begin();
	for (const auto &iter: outgoing) {
		if (ret->getPruningMultiplicity() < iter->getPruningMultiplicity())
			ret = iter;
	}
	return ret;
}

int ReadThreadingGraph::mergeDanglingHead(DanglingChainMergeHelper *danglingHeadMergeResult) {
	std::vector<CigarElement> elements = danglingHeadMergeResult->cigar->getCigarElements();
	CigarElement firstElement = elements.at(0);
	Mutect2Utils::validateArg(firstElement.getOperator() == M, "The first Cigar element must be an M");
	int indexesToMerge = bestPrefixMatch(danglingHeadMergeResult->referencePathString,
	                                     danglingHeadMergeResult->danglingPathString, firstElement.getLength());
	if (indexesToMerge <= 0) {
		return 0;
	}

	if (indexesToMerge >= danglingHeadMergeResult->referencePath.size() - 1) {
		return 0;
	}

	if (indexesToMerge >= danglingHeadMergeResult->danglingPath.size() &&
	    !extendDanglingPathAgainstReference(danglingHeadMergeResult,
	                                        indexesToMerge - danglingHeadMergeResult->danglingPath.size() + 2))
		return 0;

	addEdge(danglingHeadMergeResult->referencePath.at(indexesToMerge + 1),
	        danglingHeadMergeResult->danglingPath.at(indexesToMerge),
	        GraphObjectPool::createMultiEdge(false, 1, numPruningSamples));

	return 1;
}

int
ReadThreadingGraph::bestPrefixMatch(const std::shared_ptr<uint8_t[]> &path1, const std::shared_ptr<uint8_t[]> &path2,
                                    int maxIndex) {
	int maxMismatches = getMaxMismatches(maxIndex);
	int mismatches = 0;
	int index = 0;
	int lastGoodIndex = -1;
	uint8_t *path1_ = path1.get();
	uint8_t *path2_ = path2.get();
	while (index < maxIndex) {
		if (path1_[index] != path2_[index]) {
			if (++mismatches > maxMismatches) {
				return -1;
			}
			lastGoodIndex = index;
		}
		index++;
	}
	// if we got here then we hit the max index
	return lastGoodIndex;
}

int ReadThreadingGraph::getMaxMismatches(int lengthOfDanglingBranch) const {
	return maxMismatchesInDanglingHead > 0 ? maxMismatchesInDanglingHead : std::max(1, (lengthOfDanglingBranch /
	                                                                                    kmerSize));
}

bool ReadThreadingGraph::extendDanglingPathAgainstReference(DanglingChainMergeHelper *danglingHeadMergeResult,
                                                            int numNodesToExtend) {
	int indexOfLastDanglingNode = danglingHeadMergeResult->danglingPath.size() - 1;
	int indexOfRefNodeToUse = indexOfLastDanglingNode + numNodesToExtend;
	if (indexOfRefNodeToUse >= danglingHeadMergeResult->referencePath.size()) {
		return false;
	}

	std::shared_ptr<MultiDeBruijnVertex> danglingSource = danglingHeadMergeResult->danglingPath.at(
			indexOfLastDanglingNode);
	danglingHeadMergeResult->danglingPath.erase(
			danglingHeadMergeResult->danglingPath.begin() + indexOfLastDanglingNode);
	std::shared_ptr<uint8_t[]> danglingSourceSeq = danglingSource->getSequence();
	int danglingSourceSeqLength = danglingSource->getLength();
	std::shared_ptr<uint8_t[]> refSourceSequence = danglingHeadMergeResult->referencePath.at(
			indexOfRefNodeToUse)->getSequence();
	std::shared_ptr<uint8_t[]> sequenceToExtend(new uint8_t[numNodesToExtend + danglingSourceSeqLength]);
	uint8_t *sequenceToExtend_ = sequenceToExtend.get();
	uint8_t *refSourceSequence_ = refSourceSequence.get();
	for (int i = 0; i < numNodesToExtend; i++) {
		sequenceToExtend_[i] = refSourceSequence_[i];
	}
	memcpy(sequenceToExtend.get() + numNodesToExtend, danglingSourceSeq.get(), danglingSourceSeqLength);

	std::shared_ptr<MultiSampleEdge> sourceEdge = getHeaviestOutgoingEdge(danglingSource);
	std::shared_ptr<MultiDeBruijnVertex> prevV = getEdgeTarget(sourceEdge);
	std::shared_ptr<MultiSampleEdge> ret = removeEdge(danglingSource, prevV);
	for (int i = numNodesToExtend; i > 0; i--) {
		std::shared_ptr<uint8_t[]> tmp(new uint8_t[kmerSize]);
		memcpy(tmp.get(), sequenceToExtend.get() + i, kmerSize);
		std::shared_ptr<MultiDeBruijnVertex> newV = GraphObjectPool::createMultiVertex(tmp, kmerSize);
		addVertex(newV);
		std::shared_ptr<MultiSampleEdge> newE = addEdge(newV, prevV);
		newE->setMultiplicity(sourceEdge->getMultiplicity());
		danglingHeadMergeResult->danglingPath.emplace_back(newV);
		prevV = newV;
	}
	return true;
}

std::shared_ptr<MultiSampleEdge> ReadThreadingGraph::createEdge(const std::shared_ptr<MultiDeBruijnVertex> &,
                                                                const std::shared_ptr<MultiDeBruijnVertex> &) {
	return GraphObjectPool::createMultiEdge(false, 1, numPruningSamples);
}

std::shared_ptr<SeqGraph> ReadThreadingGraph::toSequenceGraph() {
	//buildGraphIfNecessary();
	phmap::flat_hash_set<std::shared_ptr<MultiDeBruijnVertex>> vertexSet = DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::getVertexSet();
	int reserveSize = (int) ((double) vertexSet.size() * 1.3);
	std::shared_ptr<SeqGraph> seqGraph(new SeqGraph(kmerSize));
	seqGraph->reserveSpace(reserveSize);
	phmap::flat_hash_map<std::shared_ptr<MultiDeBruijnVertex>, std::shared_ptr<SeqVertex>> vertexMap;
	vertexMap.reserve(reserveSize);
	for (auto &dv: vertexSet) {
		std::shared_ptr<SeqVertex> sv(new SeqVertex(dv->getAdditionalSequence(
				DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::isSource(dv)), dv->getAdditionalLength(
				DirectedSpecifics<MultiDeBruijnVertex, MultiSampleEdge>::isSource(dv))));
		sv->setAdditionalInfo(dv->getAdditionalInfo());
		vertexMap.insert(std::make_pair(dv, sv));
		seqGraph->addVertex(sv);
	}

	for (auto &eiter: edgeMap) {
		seqGraph->addEdge(vertexMap.at(getEdgeSource(eiter.first)), vertexMap.at(getEdgeTarget(eiter.first)),
		                  GraphObjectPool::createSeqEdge(eiter.first->getIsRef(), eiter.first->getMultiplicity()));
	}
	return seqGraph;
}

ReadThreadingGraph::ReadThreadingGraph(int kmerSize, bool debugGraphTransformations,
                                       uint8_t minBaseQualityToUseInAssembly, int numPruningSamples)
		: kmerSize(kmerSize), minBaseQualityToUseInAssembly(minBaseQualityToUseInAssembly),
		  debugGraphTransformations(debugGraphTransformations), refSource(std::make_shared<Kmer>(nullptr, 0)),
		  numPruningSamples(numPruningSamples) {
	Mutect2Utils::validateArg(kmerSize > 0, "bad minkKmerSize");
	resetToInitialState();
}

void ReadThreadingGraph::resetToInitialState() {
	pending.clear();
	nonUniqueKmers.clear();
	uniqueKmers.clear();
	alreadyBuilt = false;
}

void
ReadThreadingGraph::addSequence(std::string seqName, const std::shared_ptr<uint8_t[]> &sequence, int length,
                                int count,
                                bool isRef) {
	addSequence(std::move(seqName), ANONYMOUS_SAMPLE, sequence, 0, length, count, isRef);
}

void
ReadThreadingGraph::addSequence(std::string seqName, const std::shared_ptr<uint8_t[]> &sequence, int length,
                                bool isRef) {
	addSequence(std::move(seqName), sequence, length, 1, isRef);
}

bool ReadThreadingGraph::isLowComplexity() {
	return nonUniqueKmers.size() * 4 > uniqueKmers.size();
}

bool ReadThreadingGraph::hasCycles() {
	DFS_CycleDetect<MultiDeBruijnVertex, MultiSampleEdge> detect = DFS_CycleDetect<MultiDeBruijnVertex, MultiSampleEdge>(
			this);
	return detect.detectCycles();
}

bool ReadThreadingGraph::ifAlreadyBuilt() const {
	return alreadyBuilt;
}

void ReadThreadingGraph::printPendingInfo() {
	for (const auto &key: pendingKeys) {
		if (pending.find(key) != pending.end()) {
			std::cout << key << " " << pending[key].size() << std::endl;
			for (const auto &item: pending[key]) {
				int len = item.stop - item.start;
				char tmp[len + 1];
				memcpy(tmp, item.sequence.get() + item.start, len);
				tmp[len] = '\0';
				std::cout << item.name << std::endl << std::string(tmp) << std::endl;
			}
		}
	}
}

void ReadThreadingGraph::printGraphSize(const std::string &info) {
	if (!info.empty())
		std::cout << info << "\t";
	std::cout << "E: " << edgeMap.size() << "\t" << "V: " << vertexMapDirected.size() << "\t" << "uniqueKmers: "
	          << uniqueKmers.size() << std::endl;
}

void ReadThreadingGraph::sortPendingBySequence() {
	for (const auto &key: pendingKeys) {
		if (pending.find(key) != pending.end()) {
			std::sort(pending[key].begin(), pending[key].end(), [](SequenceForKmers &k1, SequenceForKmers &k2) -> bool {
				int len1 = k1.stop - k1.start, len2 = k2.stop - k2.start;
				if (len1 != len2)
					return len1 > len2;
				uint8_t *seq1 = k1.sequence.get() + k1.start;
				uint8_t *seq2 = k2.sequence.get() + k2.start;
				for (int i = 0; i < len1; ++i) {
					if (seq1[i] == seq2[i]) continue;
					return seq1[i] < seq2[i];
				}
				return k1.name < k2.name;
			});
		}
	}
}

