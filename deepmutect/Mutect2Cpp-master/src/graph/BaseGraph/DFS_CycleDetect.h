//
// Created by 梦想家xixi on 2021/11/5.
//

#ifndef MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H
#define MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H

#include "DirectedSpecifics.h"
#include <map>

template<class V, class E>
class DirectedSpecifics;

template<class V, class E>
class DFS_CycleDetect {
public:
	bool detectCycles() {
		phmap::flat_hash_set<std::shared_ptr<V>> vertexSet = graph->getVertexSet();
		for (const std::shared_ptr<V> &v: vertexSet) {
			if (colorMap[v] == WHITE && hasCycle(v))
				return true;
		}
		return false;
	}

	bool hasCycle(const std::shared_ptr<V> &start) {
		colorMap[start] = GREY;
		for (const std::shared_ptr<V> &v: graph->getAllTargets(start)) {
			if (colorMap[v] == GREY || (colorMap[v] == WHITE && hasCycle(v)))
				return true;
		}
		colorMap[start] = BLACK;
		return false;
	}

	explicit DFS_CycleDetect(DirectedSpecifics<V, E> *graph) : graph(graph) {
		phmap::flat_hash_set<std::shared_ptr<V>> vertexSet = graph->getVertexSet();
		for (const std::shared_ptr<V> &v: vertexSet) {
			colorMap.insert(std::make_pair(v, WHITE));
		}
	}

private:
	enum color {
		WHITE, GREY, BLACK
	};
	DirectedSpecifics<V, E> *graph;
	std::map<std::shared_ptr<V>, color> colorMap;
};


#endif //MUTECT2CPP_MASTER_DFS_CYCLEDETECT_H
