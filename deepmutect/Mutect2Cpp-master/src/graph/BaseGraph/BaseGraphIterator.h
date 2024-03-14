//
// Created by 梦想家xixi on 2021/11/16.
//

#ifndef MUTECT2CPP_MASTER_BASEGRAPHITERATOR_H
#define MUTECT2CPP_MASTER_BASEGRAPHITERATOR_H

#include <queue>
#include "DirectedSpecifics.h"
#include "parallel_hashmap/phmap.h"

template<class V, class E>
class DirectedSpecifics;

template<class T, class E>
class BaseGraphIterator {
private:
    phmap::flat_hash_set<std::shared_ptr<T>> visited;
	std::queue<std::shared_ptr<T>> toVisit;
	DirectedSpecifics<T, E> *graph;
	bool followIncomingEdges;
	bool followOutgoingEdges;

public:
	BaseGraphIterator(DirectedSpecifics<T, E> *graph, std::shared_ptr<T> start, bool followIncomingEdges,
	                  bool followOutgoingEdges) : graph(graph), followIncomingEdges(followIncomingEdges),
	                                              followOutgoingEdges(followOutgoingEdges) {
		Mutect2Utils::validateArg(graph, "graph cannot be null");
		Mutect2Utils::validateArg(start.get(), "start cannot be null");
		Mutect2Utils::validateArg(graph->containsVertex(start), "start must be in graph but it isn't");
		visited.reserve((int) (graph->getVertexSet().size() * 1.1));
		toVisit.push(start);
	}

	bool hasNext() { return !toVisit.empty(); }

	std::shared_ptr<T> next() {
		std::shared_ptr<T> v = toVisit.front();
		toVisit.pop();
		if (visited.find(v) == visited.end()) {
			visited.insert(v);
			if (followIncomingEdges) {
				for (std::shared_ptr<T> vertex: graph->vecIncomingVerticesOf(v))
					toVisit.push(vertex);
			}
			if (followOutgoingEdges) {
				for (std::shared_ptr<T> vertex: graph->vecOutgoingVerticesOf(v))
					toVisit.push(vertex);
			}
		}
		return v;
	}
};


#endif //MUTECT2CPP_MASTER_BASEGRAPHITERATOR_H
