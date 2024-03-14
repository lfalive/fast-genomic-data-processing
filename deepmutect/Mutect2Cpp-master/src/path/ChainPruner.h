//
// Created by 梦想家xixi on 2021/10/28.
//

#ifndef MUTECT2CPP_MASTER_CHAINPRUNER_H
#define MUTECT2CPP_MASTER_CHAINPRUNER_H

#include <vector>
#include "./Path.h"
#include <deque>
#include <set>


template<class V, class E>
class ChainPruner {
public:
	ChainPruner() = default;

	virtual ~ ChainPruner() = default;

	bool sortPathByStr(Path<V, E> *p1, Path<V, E> *p2) {
		int len1, len2;
		std::shared_ptr<uint8_t[]> bases1 = p1->getBases(len1);
		std::shared_ptr<uint8_t[]> bases2 = p2->getBases(len2);
		if (len1 != len2)
			return len1 > len2;
		return memcmp(bases1.get(), bases2.get(), len1) < 0;
	}

	void printAllChains(phmap::flat_hash_set<Path<V, E> *> chains) {
		std::vector<Path<V, E> *> chainsVector{};
		for (const auto &chain: chains) {
			chainsVector.push_back(chain);
		}
		printAllChains(chainsVector);
	}

	void printAllChains(std::vector<Path<V, E> *> chains) {
		if (chains.empty())
			return;
		std::cout << "chains: " << chains.size() << std::endl;
		std::vector<Path<V, E> *> sortedChains(chains.size());
		std::copy(chains.begin(), chains.end(), sortedChains.begin());
		std::sort(sortedChains.begin(), sortedChains.end(), [this](Path<V, E> *a, Path<V, E> *b) -> bool {
			return this->sortPathByStr(a, b);
		});
		int len;
		for (const auto &chain: sortedChains) {
			std::shared_ptr<uint8_t[]> bases = chain->getBases(len);
			std::cout << "edges " << chain->length() << "\t" << " bases " << len << std::endl;
			std::cout << std::string((char *) bases.get(), len) << std::endl;
		}
	}

	void pruneLowWeightChains(std::shared_ptr<DirectedSpecifics<V, E>> graph) {
		std::vector<Path<V, E> *> chains = findAllChains(graph);
		//printAllChains(chains);
		phmap::flat_hash_set<Path<V, E> *> chainsToRemoveset = chainsToRemove(chains);
		//printAllChains(chainsToRemoveset);
		//std::cout << "chains: " << chains.size() << "\t" << "remove: " << chainsToRemoveset.size() << std::endl;
		for (Path<V, E> *path: chainsToRemoveset)
			graph->removeAllEdges(path->getEdges());
		graph->removeSingletonOrphanVertices();
		for (auto chain: chains) { delete chain; }
	}

	std::vector<Path<V, E> *> findAllChains(std::shared_ptr<DirectedSpecifics<V, E>> graph) {
		std::deque<std::shared_ptr<V>> chainStarts;
		phmap::flat_hash_set<std::shared_ptr<V>> alreadySeen;
		phmap::flat_hash_set<std::shared_ptr<V>> vertexSet = graph->getVertexSet();
		for (auto &viter: vertexSet) {
			if (graph->isSource(viter)) {
				chainStarts.push_front(viter);
				alreadySeen.insert(viter);
			}
		}
		std::vector<Path<V, E> *> chains;
		while (!chainStarts.empty()) {
			std::shared_ptr<V> chainStart = chainStarts.front();
			chainStarts.pop_front();
			for (auto &eiter: graph->outgoingEdgesOf(chainStart)) {
				Path<V, E> *chain = findChain(eiter, graph);
				chains.emplace_back(chain);
				std::shared_ptr<V> chainEnd = chain->getLastVertex();
				if (alreadySeen.find(chainEnd) == alreadySeen.end()) {
					chainStarts.emplace_back(chainEnd);
					alreadySeen.insert(chainEnd);
				}
			}
		}
		return chains;
	}

private:
	Path<V, E> *findChain(std::shared_ptr<E> startEdge, std::shared_ptr<DirectedSpecifics<V, E>> graph) {
		std::vector<std::shared_ptr<E>> edges;
		edges.emplace_back(startEdge);
		std::shared_ptr<V> firstVertex = graph->getEdgeSource(startEdge);
		std::shared_ptr<V> lastVertex = graph->getEdgeTarget(startEdge);
		phmap::flat_hash_set<std::shared_ptr<E>> outEdges;
		while (true) {
			outEdges = graph->outgoingEdgesOf(lastVertex);
			if (outEdges.size() != 1 || graph->inDegreeOf(lastVertex) > 1 || lastVertex == firstVertex) {
				break;
			}
			std::shared_ptr<E> nextEdge = *outEdges.begin();
			edges.emplace_back(nextEdge);
			lastVertex = graph->getEdgeTarget(nextEdge);
		}
		return new Path<V, E>(edges, lastVertex, graph);
	}

protected:
	virtual phmap::flat_hash_set<Path<V, E> *> chainsToRemove(std::vector<Path<V, E> *> chains) = 0;
};


#endif //MUTECT2CPP_MASTER_CHAINPRUNER_H
