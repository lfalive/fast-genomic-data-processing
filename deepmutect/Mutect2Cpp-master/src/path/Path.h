//
// Created by 梦想家xixi on 2021/10/27.
//

#ifndef MUTECT2CPP_MASTER_PATH_H
#define MUTECT2CPP_MASTER_PATH_H

#include <cstring>
#include <vector>
#include "graph/BaseGraph/DirectedSpecifics.h"
#include "Mutect2Utils.h"
#include <iostream>

template<class T, class E>
class Path {
private:
	// the last vertex seen in the path
	std::shared_ptr<T> lastVertex;

	// the list of edges comprising the path
	std::vector<std::shared_ptr<E>> edgesInOrder;

	// the graph from which this path originated
	std::shared_ptr<DirectedSpecifics<T, E>> graph;

    double logOdds;

public:
	/**
	 * Create a new Path containing no edges and starting at initialVertex
	 * @param initialVertex the starting vertex of the path
	 * @param graph the graph this path will follow through
	 */
	Path(std::shared_ptr<T> initialVertex, std::shared_ptr<DirectedSpecifics<T, E>> graph) : lastVertex(initialVertex),
	                                                                                         graph(graph) {
		Mutect2Utils::validateArg(initialVertex.get(), "initialVertex cannot be null");
		Mutect2Utils::validateArg(graph->containsVertex(initialVertex), "Vertex must be part of graph.");
	}

	/**
	 * Constructor that does not check arguments' validity i.e. doesn't check that edges are in order
	 */
	Path(std::vector<std::shared_ptr<E>> edgesInOrder, std::shared_ptr<T> lastVertex,
	     std::shared_ptr<DirectedSpecifics<T, E>> graph) : lastVertex(lastVertex), graph(graph),
	                                                       edgesInOrder(edgesInOrder) {}

	/**
	 * Create a new Path extending p with edge
	 *
	 * @param p the path to extend.
	 * @param edge the edge to extend path with.
	 *
	 * @throws IllegalArgumentException if {@code p} or {@code edge} are {@code null}, or {@code edge} is
	 * not part of {@code p}'s graph, or {@code edge} does not have as a source the last vertex in {@code p}.
	 */
	Path(Path<T, E> p, std::shared_ptr<E> edge) : graph(p.graph), lastVertex(p.graph->getEdgeTarget(edge)) {
		Mutect2Utils::validateArg(edge.get(), "Edge cannot be null");
		Mutect2Utils::validateArg(p.graph->getEdgeSource(edge) == p.lastVertex,
		                          "Edges added to path must be contiguous.");
		for (typename std::vector<std::shared_ptr<E>>::iterator iter = p.edgesInOrder.begin();
		     iter != p.edgesInOrder.end(); iter++) {
			edgesInOrder.template emplace_back(*iter);
		}
		edgesInOrder.template emplace_back(edge);
	}

	int length() const { return edgesInOrder.size(); }

	Path(std::shared_ptr<E> edge, Path<T, E> p) : graph(p.graph), lastVertex(p.lastVertex) {
		Mutect2Utils::validateArg(edge, "Edge cannot be null");
		Mutect2Utils::validateArg(p.graph.containsEdge(edge), "Graph must contain edge ");
		Mutect2Utils::validateArg(p.graph.getEdgeTarget(edge) == p.getFirstVertex(),
		                          "Edges added to path must be contiguous");
		for (typename std::vector<std::shared_ptr<E>>::iterator iter = p.edgesInOrder.begin();
		     iter != p.edgesInOrder.end(); iter++) {
			edgesInOrder.insert(*iter);
		}
		edgesInOrder.insert(edge);
	}

	bool containsVertex(std::shared_ptr<T> v) {
		Mutect2Utils::validateArg(v, "Vertex cannot be null");
		std::vector<std::shared_ptr<T>> res;
		res.emplace_back(getFirstVertex());
		for (typename std::vector<std::shared_ptr<E>>::iterator iter = edgesInOrder.begin();
		     iter != edgesInOrder.end(); iter++) {
			res.emplace_back(graph->getEdgeTarget(*iter));
		}
		return std::find(res.begin(), res.end(), v) != res.end();
	}

	std::shared_ptr<DirectedSpecifics<T, E>> getGraph() { return graph; }

	std::vector<std::shared_ptr<E>> &getEdges() { return edgesInOrder; }

	std::shared_ptr<E> getLastEdge() { return edgesInOrder[edgesInOrder.size() - 1]; }

	std::vector<std::shared_ptr<T>> getVertices() {
		std::vector<std::shared_ptr<T>> res;
		res.emplace_back(getFirstVertex());
		for (typename std::vector<std::shared_ptr<E>>::iterator iter = edgesInOrder.begin();
		     iter != edgesInOrder.end(); iter++) {
			res.emplace_back(graph->getEdgeTarget(*iter));
		}
		return res;
	}

	std::shared_ptr<T> getFirstVertex() {
		if (edgesInOrder.empty()) {
			return lastVertex;
		} else {
			return getGraph()->getEdgeSource(edgesInOrder[0]);
		}
	}

	std::shared_ptr<uint8_t[]> getBases(int &reslength) {
		if (getEdges().empty()){
			reslength = graph->getAdditionalSequenceLength(lastVertex, true);
			return graph->getAdditionalSequence(lastVertex, true);
		}

		std::shared_ptr<T> source = graph->getEdgeSource(edgesInOrder[0]);
		std::shared_ptr<uint8_t[]> bases = graph->getAdditionalSequence(source, true);
		int basesLength = graph->getAdditionalSequenceLength(source, true);

		std::vector<uint8_t> baseVec(bases.get(), bases.get() + basesLength);

		for (int i = 0; i < edgesInOrder.size(); i++) {
			std::shared_ptr<T> target = graph->getEdgeTarget(edgesInOrder[i]);
			bases = graph->getAdditionalSequence(target, false);
			basesLength = graph->getAdditionalSequenceLength(target, false);

			for (int k = 0; k < basesLength; ++k) {
				baseVec.push_back(bases.get()[k]);
			}
		}
		reslength = (int) baseVec.size();
		std::shared_ptr<uint8_t[]> ret_bases(new uint8_t[reslength]);
		memcpy(ret_bases.get(), &baseVec[0], reslength);

		return ret_bases;
	}

	std::shared_ptr<T> getLastVertex() { return lastVertex; }

    double getLogOdds() {return logOdds;}

    void setLogOdds(double score) {
        this->logOdds = score;
    }
};

#endif //MUTECT2CPP_MASTER_PATH_H
