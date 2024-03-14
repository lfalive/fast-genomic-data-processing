//
// Created by 梦想家xixi on 2021/11/19.
//

#include "CommonSuffixSplitter.h"
#include <memory>
#include "parallel_hashmap/phmap.h"
#include "GraphUtils.h"
#include "GraphObjectPool.h"

bool CommonSuffixSplitter::split(SeqGraph *graph, std::shared_ptr<SeqVertex> v) {
	Mutect2Utils::validateArg(graph, "graph cannot be null");
	Mutect2Utils::validateArg(v.get(), "v cannot be null");
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &allVertex = graph->getVertexSet();
	Mutect2Utils::validateArg(allVertex.find(v) != allVertex.end(), "graph doesn't contain vertex v ");
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> toSplit = graph->incomingVerticesOf(v);
	std::shared_ptr<SeqVertex> suffixVTemplate = commonSuffix(graph, v, toSplit);
	if (suffixVTemplate == nullptr) {
		return false;
	}
	std::list<std::shared_ptr<BaseEdge>> edgesToRemove;
#ifdef SORT_MODE
	for (const std::shared_ptr<SeqVertex> &mid: graph->sortedVerticesOf(toSplit)) {
#else
	for (const std::shared_ptr<SeqVertex> &mid: toSplit) {
#endif
		std::shared_ptr<SeqVertex> suffixV = GraphObjectPool::createSeqVertex(suffixVTemplate->getSequence(),
		                                                                 suffixVTemplate->getLength());
		graph->addVertex(suffixV);
		std::shared_ptr<SeqVertex> prefixV = mid->withoutSuffix(suffixV->getSequence(), suffixV->getLength());
		std::shared_ptr<BaseEdge> out = graph->outgoingEdgeOf(mid);
		std::shared_ptr<SeqVertex> incomingTarget;
		if (prefixV == nullptr) {
			incomingTarget = suffixV;
		} else {
			incomingTarget = prefixV;
			graph->addVertex(prefixV);
			graph->addEdge(prefixV, suffixV, GraphObjectPool::createSeqEdge(out->getIsRef(), 1));
			edgesToRemove.emplace_back(out);
		}
		graph->addEdge(suffixV, graph->getEdgeTarget(out),
		               GraphObjectPool::createSeqEdge(out->getIsRef(), out->getMultiplicity()));
		for (const std::shared_ptr<BaseEdge> &in: graph->incomingEdgesOf(mid)) {
			graph->addEdge(graph->getEdgeSource(in), incomingTarget,
			               GraphObjectPool::createSeqEdge(in->getIsRef(), in->getMultiplicity()));
			edgesToRemove.emplace_back(in);
		}
	}
	graph->removeAllVertices(toSplit);
	graph->removeAllEdges(edgesToRemove);
	return true;
}

std::shared_ptr<SeqVertex> CommonSuffixSplitter::commonSuffix(SeqGraph *graph, std::shared_ptr<SeqVertex> v,
                                                              const phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &toSplit) {
	if (toSplit.size() < 2) {
		return nullptr;
	} else if (!safeToSplit(graph, v, toSplit)) {
		return nullptr;
	}
	std::shared_ptr<SeqVertex> suffixVTemplate = commonSuffix(toSplit);
	if (suffixVTemplate->isEmpty()) {
		return nullptr;
	} else if (wouldEliminateRefSource(graph, suffixVTemplate, toSplit)) {
		return nullptr;
	} else if (allVerticesAreTheCommonSuffix(suffixVTemplate, toSplit)) {
		return nullptr;
	} else {
		return suffixVTemplate;
	}
}

bool CommonSuffixSplitter::safeToSplit(SeqGraph *graph, const std::shared_ptr<SeqVertex> &bot,
                                       const phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &toMerge) {
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> outgoingVertices = graph->outgoingVerticesOf(bot);
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> outgoingOfBot;
	for (const std::shared_ptr<SeqVertex> &v: outgoingVertices) {
		outgoingOfBot.insert(v);
	}
	for (const std::shared_ptr<SeqVertex> &m: toMerge) {
		phmap::flat_hash_set<std::shared_ptr<BaseEdge>> outs = graph->outgoingEdgesOf(m);
		phmap::flat_hash_set<std::shared_ptr<SeqVertex>> tmp = graph->outgoingVerticesOf(m);
		if (m == bot || outs.size() != 1 || tmp.find(bot) == tmp.end()) {
			return false;
		}
		if (outgoingOfBot.find(m) != outgoingOfBot.end()) {
			return false;
		}
	}
	return true;
}

std::shared_ptr<SeqVertex>
CommonSuffixSplitter::commonSuffix(const phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &middleVertices) {
	std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> kmers = GraphUtils::getKmers(middleVertices);
	int min = GraphUtils::minKmerLength(kmers);
	int suffixLen = GraphUtils::commonMaximumSuffixLength(kmers, min);
	std::shared_ptr<uint8_t[]> kmer = kmers.begin()->first;
	int kmerLength = kmers.begin()->second;
	int suffixLength;
	std::shared_ptr<uint8_t[]> suffix = Mutect2Utils::copyOfRange(kmer, kmerLength, kmerLength - suffixLen, kmerLength,
	                                                              suffixLength);
	return GraphObjectPool::createSeqVertex(suffix, suffixLength);
}

bool
CommonSuffixSplitter::wouldEliminateRefSource(SeqGraph *graph, const std::shared_ptr<SeqVertex> &commonSuffix,
                                              const phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &toSplits) {
	for (const std::shared_ptr<SeqVertex> &toSplit: toSplits) {
		if (graph->isRefSource(toSplit)) {
			return toSplit->getLength() == commonSuffix->getLength();
		}
	}
	return false;
}

bool CommonSuffixSplitter::allVerticesAreTheCommonSuffix(const std::shared_ptr<SeqVertex> &commonSuffix,
                                                         const phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &toSplits) {
	for (const std::shared_ptr<SeqVertex> &toSplit: toSplits) {
		if (toSplit->getLength() != commonSuffix->getLength()) {
			return false;
		}
	}
	return true;
}
