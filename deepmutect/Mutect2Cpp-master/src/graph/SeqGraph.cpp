//
// Created by 梦想家xixi on 2021/11/15.
//

#include "SeqGraph.h"
#include <list>
#include <climits>
#include <cstring>
#include <memory>
#include <algorithm>
#include "utils/MergeCommonSuffices.h"
#include "utils/MergeDiamonds.h"
#include "utils/MergeTails.h"
#include "utils/SplitCommonSuffices.h"
#include "utils/GraphUtils.h"
#include "graph/utils/GraphObjectPool.h"

bool SeqGraph::zipLinearChains() {
	std::vector<std::shared_ptr<SeqVertex>> zipStarts;
#ifdef SORT_MODE
	std::vector<std::shared_ptr<SeqVertex>> vertexSet = getSortedVertexList();
#else
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> vertexSet = getVertexSet();
#endif
	for (auto &source: vertexSet) {
		if (isLinearChainStart(source)) {
			zipStarts.emplace_back(source);
		}
	}

	if (zipStarts.empty())
		return false;

	std::list<std::shared_ptr<SeqVertex>> linearChain;
	bool mergedOne = false;
	for (auto &zipStart: zipStarts) {
		linearChain = traceLinearChain(zipStart);
		mergedOne |= mergeLinearChain(linearChain);
	}
	return mergedOne;
}

bool SeqGraph::isLinearChainStart(const std::shared_ptr<SeqVertex> &source) {
	return outDegreeOf(source) == 1
	       && (inDegreeOf(source) != 1 || outDegreeOf(*(incomingVerticesOf(source).begin())) > 1);
}

std::list<std::shared_ptr<SeqVertex>> SeqGraph::traceLinearChain(const std::shared_ptr<SeqVertex> &zipStart) {
	std::list<std::shared_ptr<SeqVertex>> linearChain;
	linearChain.emplace_back(zipStart);

	bool lastIsRef = isReferenceNode(zipStart);
	std::shared_ptr<SeqVertex> last = zipStart;
	while (true) {
		if (outDegreeOf(last) != 1)
			break;

		std::shared_ptr<SeqVertex> target = getEdgeTarget(outgoingEdgeOf(last));
		if (inDegreeOf(target) != 1 || last == target)
			break;

		bool targetIsRef = isReferenceNode(target);
		if (lastIsRef != targetIsRef)
			break;
		linearChain.emplace_back(target);
		last = target;
		lastIsRef = targetIsRef;
	}
	return linearChain;
}

bool SeqGraph::mergeLinearChain(std::list<std::shared_ptr<SeqVertex>> &linearChain) {
	if (linearChain.empty())
		throw std::invalid_argument("BUG: cannot have linear chain with 0 elements");

	std::shared_ptr<SeqVertex> first = linearChain.front();
	std::shared_ptr<SeqVertex> last = linearChain.back();

	if (first == last)
		return false;

	std::shared_ptr<SeqVertex> addedVertex = mergeLinearChainVertices(linearChain);
	addVertex(addedVertex);

	for (auto &edge: outgoingEdgesOf(last)) {
		addEdge(addedVertex, getEdgeTarget(edge),
		        GraphObjectPool::createSeqEdge(edge->getIsRef(), edge->getMultiplicity()));
	}

	for (auto &edge: incomingEdgesOf(first)) {
		addEdge(getEdgeSource(edge), addedVertex,
		        GraphObjectPool::createSeqEdge(edge->getIsRef(), edge->getMultiplicity()));
	}
	removeAllVertices(linearChain);
	return true;
}


std::shared_ptr<SeqVertex> SeqGraph::mergeLinearChainVertices(std::list<std::shared_ptr<SeqVertex>> &vertices) {
	int length = 500;
	int start = 0;
	std::shared_ptr<uint8_t[]> tmp(new uint8_t[length]);
	for (auto &v: vertices) {
		int seqLength = v->getLength();
		while (start + seqLength >= length) {
			length <<= 1;
			std::shared_ptr<uint8_t[]> newtmp(new uint8_t[length]);
			memcpy(newtmp.get(), tmp.get(), start);
			tmp = newtmp;
		}
		memcpy(tmp.get() + start, v->getSequence().get(), seqLength);
		start += seqLength;
	}
	return GraphObjectPool::createSeqVertex(tmp, start);
}

void SeqGraph::simplifyGraph() {
	simplifyGraph(INT_MAX);
}

void SeqGraph::simplifyGraph(int maxCycles) {
	zipLinearChains();
	std::shared_ptr<SeqGraph> prevGraph = nullptr;
	for (int i = 0; i < maxCycles; i++) {
		if (i > MAX_REASONABLE_SIMPLIFICATION_CYCLES) {
			throw std::invalid_argument("Infinite loop detected in simplification routines for kmer graph");
		}
		if (!simplifyGraphOnce(i))
			break;
		if (i > 5) {
			if (prevGraph != nullptr && GraphUtils::graphEquals(prevGraph.get(), this))
				break;
			prevGraph = std::shared_ptr<SeqGraph>(clone());
		}
	}
}

bool SeqGraph::simplifyGraphOnce(int iteration) {
	bool didSomeWork = MergeDiamonds(this).transformUntilComplete();
	//printGraphSize(std::to_string(iteration) + " MergeDiamonds");
	//outputDotFile("CPP_MergeDiamonds.dot");

	didSomeWork |= MergeTails(this).transformUntilComplete();
	//printGraphSize(std::to_string(iteration) + " MergeTails");
	//outputDotFile("CPP_MergeTails.dot");

	didSomeWork |= SplitCommonSuffices(this).transformUntilComplete();
	//printGraphSize(std::to_string(iteration) + " SplitCommonSuffices");
	//outputDotFile("CPP_SplitCommonSuffices.dot");

	didSomeWork |= MergeCommonSuffices(this).transformUntilComplete();
	//printGraphSize(std::to_string(iteration) + " MergeCommonSuffices");
	//outputDotFile("CPP_MergeCommonSuffices.dot");

	didSomeWork |= zipLinearChains();
	//printGraphSize(std::to_string(iteration) + " zipLinearChains");
	//outputDotFile("CPP_zipLinearChains.dot");
	return didSomeWork;
}

SeqGraph::SeqGraph(SeqGraph &seqGraph) : kmerSize(seqGraph.kmerSize),
                                         DirectedSpecifics<SeqVertex, BaseEdge>(seqGraph.getVertexSet(),
                                                                                seqGraph.getEdgeSet()) {
	vertexMapDirected = seqGraph.vertexMapDirected;
	edgeMap = seqGraph.edgeMap;
}

SeqGraph *SeqGraph::clone() {
	return new SeqGraph(*this);
}

void SeqGraph::printGraphSize(const std::string &info) {
	if (!info.empty())
		std::cout << info << "\t";
	std::cout << "E: " << edgeMap.size() << "\t" << "V: " << vertexMapDirected.size() << std::endl;
}

