//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_GRAPHUTILS_H
#define MUTECT2CPP_MASTER_GRAPHUTILS_H

#include <list>
#include <set>
#include "SeqVertex.h"
#include "BaseGraph/DirectedSpecifics.h"

class GraphUtils {
public:
	static int minKmerLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> &kmers);

	static int commonMaximumPrefixLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> &kmers);

	static int
	commonMaximumSuffixLength(std::list<std::pair<std::shared_ptr<uint8_t[]>, int>> &listOfBytes, int minLength);

	static std::list<std::pair<std::shared_ptr<uint8_t[]>, int>>
	getKmers(const std::vector<std::shared_ptr<SeqVertex>> &vertices);

	static std::list<std::pair<std::shared_ptr<uint8_t[]>, int>>
	getKmers(const phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &vertices);

	template<class T, class E>
	static bool graphEquals(DirectedSpecifics<T, E> *g1, DirectedSpecifics<T, E> *g2) {
		Mutect2Utils::validateArg(g1, "g1");
		Mutect2Utils::validateArg(g2, "g2");
		phmap::flat_hash_set<std::shared_ptr<T>> vertices1 = g1->getVertexSet();
		phmap::flat_hash_set<std::shared_ptr<T>> vertices2 = g2->getVertexSet();
		phmap::flat_hash_set<std::shared_ptr<E>> edges1 = g1->getEdgeSet();
		phmap::flat_hash_set<std::shared_ptr<E>> edges2 = g2->getEdgeSet();
		if (vertices1.size() != vertices2.size() || edges1.size() != edges2.size())
			return false;

		// Check that for every vertex in g1 there is a vertex in g2 with an equal getSequenceString
		bool all_vertices_match = std::all_of(vertices1.begin(), vertices1.end(), [&](const std::shared_ptr<T> &v1) {
			return std::any_of(vertices2.begin(), vertices2.end(), [&](const std::shared_ptr<T> &v2) {
				return (*v1 == *v2);
			});
		});
		if (!all_vertices_match)
			return false;

		// Check that for every edge in g1 there is an equal edge in g2
		bool all_edges_match_g1 = std::all_of(edges1.begin(), edges1.end(), [&](const std::shared_ptr<E> &edge1) {
			return std::any_of(edges2.begin(), edges2.end(), [&](const std::shared_ptr<E> &edge2) {
				return (*g1->getEdgeSource(edge1) == *g2->getEdgeSource(edge2)) &&
				       (*g1->getEdgeTarget(edge1) == *g2->getEdgeTarget(edge2));
			});
		});
		if (!all_edges_match_g1)
			return false;

		// Check that for every edge in g2 there is an equal edge in g1
		bool all_edges_match_g2 = std::all_of(edges2.begin(), edges2.end(), [&](const std::shared_ptr<E> &edge2) {
			return std::any_of(edges1.begin(), edges1.end(), [&](const std::shared_ptr<E> &edge1) {
				return (*g1->getEdgeSource(edge1) == *g2->getEdgeSource(edge2)) &&
				       (*g1->getEdgeTarget(edge1) == *g2->getEdgeTarget(edge2));
			});
		});
		return all_edges_match_g2;
	};
};


#endif //MUTECT2CPP_MASTER_GRAPHUTILS_H
