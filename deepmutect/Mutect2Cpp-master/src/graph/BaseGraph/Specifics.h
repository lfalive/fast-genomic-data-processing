//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_SPECIFICS_H
#define MUTECT2CPP_MASTER_SPECIFICS_H

#include <memory>
#include <stdexcept>
#include <vector>
#include "parallel_hashmap/phmap_fwd_decl.h"

template<class V, class E>
class Specifics {
public:
	virtual void addVertex(const std::shared_ptr<V> &vertex) = 0;

	virtual phmap::flat_hash_set<std::shared_ptr<V>> &getVertexSet() = 0;

	virtual phmap::flat_hash_set<std::shared_ptr<E>>
	getAllEdges(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) = 0;

	virtual std::shared_ptr<E>
	getEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) = 0;

	virtual phmap::flat_hash_set<std::shared_ptr<E>> edgesof(const std::shared_ptr<V> &vertex) = 0;

	virtual int inDegreeOf(const std::shared_ptr<V> &vertex) = 0;

	virtual phmap::flat_hash_set<std::shared_ptr<E>> &incomingEdgesOf(const std::shared_ptr<V> &vertex) = 0;

	virtual int outDegreeOf(const std::shared_ptr<V> &vertex) = 0;

	virtual phmap::flat_hash_set<std::shared_ptr<E>> &outgoingEdgesOf(const std::shared_ptr<V> &vertex) = 0;

	virtual void removeEdgeFromTouchingVertices(const std::shared_ptr<E> &e) = 0;

	virtual bool isSource(const std::shared_ptr<V> &v) = 0;

	virtual bool isSink(const std::shared_ptr<V> &v) = 0;

	virtual bool isSingletonOrphan(const std::shared_ptr<V> &v) {
		if (!v.get())
			throw std::invalid_argument("v can not be null.");
		return inDegreeOf(v) == 0 && outDegreeOf(v) == 0 && !isRefSource(v);
	}

	virtual bool isRefSource(const std::shared_ptr<V> &v) {
		if (!v.get())
			throw std::invalid_argument("Attempting to pull sequence from a null vertex.");
		for (auto &e: incomingEdgesOf(v)) {
			if (e->getIsRef())
				return false;
		}
		for (auto &e: outgoingEdgesOf(v)) {
			if (e->getIsRef())
				return true;
		}

		return getVertexSet().size() == 1;
	}

	virtual void removeSingletonOrphanVertices() {
		std::vector<std::shared_ptr<V>> toRemove;
		phmap::flat_hash_set<std::shared_ptr<V>> &allvertex = getVertexSet();
		typename phmap::flat_hash_set<std::shared_ptr<V>>::iterator viter;
		for (viter = allvertex.begin(); viter != allvertex.end(); viter++) {
			if (isSingletonOrphan(*viter)) {
				toRemove.template emplace_back(*viter);
			}
		}
		removeAllVertices(toRemove);
	}

	virtual bool removeAllVertices(const std::vector<std::shared_ptr<V>> &vertices) = 0;

	virtual bool removeAllEdges(const std::vector<std::shared_ptr<E>> &edges) = 0;
};


#endif //MUTECT2CPP_MASTER_SPECIFICS_H
