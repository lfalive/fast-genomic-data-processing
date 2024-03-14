//
// Created by 梦想家xixi on 2021/11/20.
//

#include <cstring>
#include <memory>
#include "SharedSequenceMerger.h"
#include "graph/utils/GraphObjectPool.h"

bool SharedSequenceMerger::canMerge(SeqGraph *graph, const std::shared_ptr<SeqVertex> &v,
                                    phmap::flat_hash_set<std::shared_ptr<SeqVertex>> incomingVertices) {
	if (incomingVertices.empty()) {
		return false;
	}

	std::shared_ptr<SeqVertex> first = *incomingVertices.begin();
	for (const std::shared_ptr<SeqVertex> &prev: incomingVertices) {
		if (!prev->seqEquals(first)) {
			return false;
		}
		std::vector<std::shared_ptr<SeqVertex>> prevOuts = graph->vecOutgoingVerticesOf(prev);
		if (prevOuts.size() != 1) {
			return false;
		}
		if (*prevOuts.begin() != v) {
			return false;
		}
		if (graph->inDegreeOf(prev) == 0) {
			return false;
		}
	}
	return true;
}

bool SharedSequenceMerger::merge(SeqGraph *graph, const std::shared_ptr<SeqVertex> &v) {
	Mutect2Utils::validateArg(graph, "graph cannot be null");
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &allVertex = graph->getVertexSet();
	Mutect2Utils::validateArg(allVertex.find(v) != allVertex.end(), "graph doesn't contain vertex");
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> prevs = graph->incomingVerticesOf(v);
	if (!canMerge(graph, v, prevs)) {
		return false;
	} else {
		std::list<std::shared_ptr<BaseEdge>> edgesToRemove;
		std::shared_ptr<uint8_t[]> prevSeq = (*prevs.begin())->getSequence();
		int prevSeqLength = (*prevs.begin())->getLength();
		std::shared_ptr<uint8_t[]> vSeq = v->getSequence();
		int vSeqLength = v->getLength();
		int tmpLength = prevSeqLength + vSeqLength;
		std::shared_ptr<uint8_t[]> tmp(new uint8_t[tmpLength]);
		memcpy(tmp.get(), prevSeq.get(), prevSeqLength);
		memcpy(tmp.get() + prevSeqLength, vSeq.get(), vSeqLength);
		std::shared_ptr<SeqVertex> newV(new SeqVertex(tmp, tmpLength));
		graph->addVertex(newV);
#ifdef SORT_MODE
		for (const std::shared_ptr<SeqVertex> &prev: graph->sortedVerticesOf(prevs)) {
#else
		for (const std::shared_ptr<SeqVertex> &prev: prevs) {
#endif
			for (const std::shared_ptr<BaseEdge> &prevIn: graph->incomingEdgesOf(prev)) {
				graph->addEdge(graph->getEdgeSource(prevIn), newV,
				               GraphObjectPool::createSeqEdge(prevIn->getIsRef(), prevIn->getMultiplicity()));
				edgesToRemove.emplace_back(prevIn);
			}
		}
		for (const std::shared_ptr<BaseEdge> &e: graph->outgoingEdgesOf(v)) {
			graph->addEdge(newV, graph->getEdgeTarget(e),
			               GraphObjectPool::createSeqEdge(e->getIsRef(), e->getMultiplicity()));
		}
		graph->removeAllVertices(prevs);
		graph->removeVertex(v);
		graph->removeAllEdges(edgesToRemove);
		return true;
	}
}
