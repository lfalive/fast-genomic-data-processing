//
// Created by 梦想家xixi on 2021/10/21.
//

#ifndef MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
#define MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H

#include "Specifics.h"
#include "Mutect2Utils.h"
#include "DirectedEdgeContainer.h"
#include "BaseGraphIterator.h"
#include <stdexcept>
#include <cstring>
#include <deque>
#include <list>
#include <map>
#include <fstream>
#include <algorithm>
#include "DFS_CycleDetect.h"

static const std::string LOOPS_NOT_ALLOWED = "loops not allowed";

template<class T, class E>
class BaseGraphIterator;


template<class V>
class IntrusiveEdge {
private:
	std::shared_ptr<V> source;
	std::shared_ptr<V> target;

public:
	IntrusiveEdge(std::shared_ptr<V> source, std::shared_ptr<V> target) : source(source), target(target) {};

	std::shared_ptr<V> getSource() const { return source; }

	std::shared_ptr<V> getTarget() const { return target; }
};

template<class V, class E>
class DirectedSpecifics : public Specifics<V, E> {
private:
	DirectedEdgeContainer<E> &getEdgeContainer(const std::shared_ptr<V> &vertex) {
		auto miter = vertexMapDirected.find(vertex);
		if (miter == vertexMapDirected.end()) {
			throw std::invalid_argument("no such vertex in graph");
		}
		return miter->second;
	}

	bool allowingLoops = true;
	bool allowingMultipleEdges = false;
	phmap::flat_hash_set<std::shared_ptr<V>> VertexSet;
	phmap::flat_hash_set<std::shared_ptr<E>> EdgeSet;

public:
	phmap::flat_hash_map<std::shared_ptr<V>, DirectedEdgeContainer<E>> vertexMapDirected;

	phmap::flat_hash_map<std::shared_ptr<E>, IntrusiveEdge<V>> edgeMap;

	DirectedSpecifics() = default;

	DirectedSpecifics(const phmap::flat_hash_set<std::shared_ptr<V>> &vertexSet,
	                  const phmap::flat_hash_set<std::shared_ptr<E>> &edgeSet) : VertexSet(vertexSet), EdgeSet(edgeSet) {}

	~DirectedSpecifics() = default;

	int getEdgesNum() {
		return EdgeSet.size();
	}

	int getVertexNum() {
		return VertexSet.size();
	}

	void addVertex(const std::shared_ptr<V> &v) {
		if (v == nullptr)
			throw std::invalid_argument("Null is not allowed here.");
		if (containsVertex(v))
			return;
		vertexMapDirected.insert(std::make_pair(v, DirectedEdgeContainer<E>()));
		VertexSet.insert(v);
	}

	std::vector<std::shared_ptr<V>> sortedVerticesOf(phmap::flat_hash_set<std::shared_ptr<V>> vertices) {
		std::vector<std::shared_ptr<V>> ret;
		ret.reserve(vertices.size());
		for (auto &v: vertices) {
			ret.emplace_back(v);
		}
		std::sort(ret.begin(), ret.end(), [this](std::shared_ptr<V> v1, std::shared_ptr<V> v2) -> bool {
			int len1 = v1->getLength();
			int len2 = v2->getLength();
			if (len1 != len2)
				return len1 > len2;
			uint8_t *seq1 = v1->getSequence().get();
			uint8_t *seq2 = v2->getSequence().get();
			for (int i = 0; i < len1; ++i) {
				if (seq1[i] == seq2[i]) continue;
				return seq1[i] < seq2[i];
			}

			phmap::flat_hash_set<std::shared_ptr<E>> in1 = incomingEdgesOf(v1);
			phmap::flat_hash_set<std::shared_ptr<E>> in2 = incomingEdgesOf(v2);
			if (in1.size() != in2.size())
				return in1.size() > in2.size();

			phmap::flat_hash_set<std::shared_ptr<E>> out1 = outgoingEdgesOf(v1);
			phmap::flat_hash_set<std::shared_ptr<E>> out2 = outgoingEdgesOf(v2);
			if (out1.size() != out2.size())
				return out1.size() > out2.size();

			int val1 = 0, val2 = 0;
			for (auto &e: in1)
				val1 += getEdgeSource(e)->getLength() * e->getMultiplicity();
			for (auto &e: in2)
				val2 += getEdgeSource(e)->getLength() * e->getMultiplicity();
			if (val1 != val2)
				return val1 < val2;

			val1 = 0, val2 = 0;
			for (auto &e: out1)
				val1 += getEdgeTarget(e)->getLength() * e->getMultiplicity();
			for (auto &e: out2)
				val2 += getEdgeTarget(e)->getLength() * e->getMultiplicity();
			if (val1 != val2)
				return val1 < val2;

			val1 = 0, val2 = 0;
			for (auto &e: in1)
				val1 += getEdgeSource(e)->getSequence().get()[0] * e->getMultiplicity();
			for (auto &e: in2)
				val2 += getEdgeSource(e)->getSequence().get()[0] * e->getMultiplicity();
			if (val1 != val2)
				return val1 < val2;

			val1 = 0, val2 = 0;
			for (auto &e: out1)
				val1 += getEdgeTarget(e)->getSequence().get()[0] * e->getMultiplicity();
			for (auto &e: out2)
				val2 += getEdgeTarget(e)->getSequence().get()[0] * e->getMultiplicity();
			return val1 < val2;
		});
		return ret;
	}

	std::vector<std::shared_ptr<V>> getSortedVertexList() {
		return sortedVerticesOf(VertexSet);
	}

	phmap::flat_hash_set<std::shared_ptr<V>> &getVertexSet() {
		return VertexSet;
	}

	phmap::flat_hash_set<std::shared_ptr<E>>
	getAllEdges(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		phmap::flat_hash_set<std::shared_ptr<E>> edges;
		if (VertexSet.find(sourceVertex) != VertexSet.end() && VertexSet.find(targetVertex) != VertexSet.end()) {
			for (auto &edge: getEdgeContainer(sourceVertex).outgoing) {
				if (getEdgeTarget(edge) == targetVertex)
					edges.insert(edge);
			}
		}
		return edges;
	}

	phmap::flat_hash_set<std::shared_ptr<V>> getAllTargets(const std::shared_ptr<V> &sourceVertex) {
		phmap::flat_hash_set<std::shared_ptr<V>> res;
		if (VertexSet.find(sourceVertex) != VertexSet.end()) {
			for (auto &vertex: getEdgeContainer(sourceVertex).outgoing) {
				res.insert(getEdgeTarget(vertex));
			}
		}
		return res;
	}

	std::shared_ptr<V> getEdgeTarget(const std::shared_ptr<E> &e) {
		return edgeMap.find(e)->second.getTarget();
	}

	std::shared_ptr<E> getEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		if (VertexSet.find(sourceVertex) != VertexSet.end() && VertexSet.find(targetVertex) != VertexSet.end()) {
			for (auto &edge: getEdgeContainer(sourceVertex).outgoing) {
				if (getEdgeTarget(edge) == targetVertex)
					return edge;
			}
		}
		return nullptr;
	}

	std::shared_ptr<V> getEdgeSource(std::shared_ptr<E> e) {
		return edgeMap.find(e)->second.getSource();
	}

	phmap::flat_hash_set<std::shared_ptr<E>> edgesof(const std::shared_ptr<V> &vertex) {
		phmap::flat_hash_set<std::shared_ptr<E>> res = getEdgeContainer(vertex).incoming;
		const phmap::flat_hash_set<std::shared_ptr<E>> &outgoing = getEdgeContainer(vertex).outgoing;
		res.reserve(res.size() + outgoing.size());
		for (const std::shared_ptr<E> &e: outgoing) {
			res.insert(e);
		}

		if (allowingLoops) {
			phmap::flat_hash_set<std::shared_ptr<E>> loops = getAllEdges(vertex, vertex);
			for (auto &v: res) {
				if (loops.find(v) != loops.end()) {
					loops.erase(v);
					res.erase(v);
				}
			}
		}
		return res;
	}

	void outputDotFile(const std::string &fileName) {
		std::ofstream outfile1(fileName);
		outfile1 << "digraph G{" << std::endl;
		for (auto &v: getSortedVertexList()) {
			std::string s(reinterpret_cast<const char *>(v->getSequence().get()), v->getLength());
			outfile1 << "    " << s;
			if (isReferenceNode(v))
				outfile1 << "[color=Red]";
			outfile1 << ";" << std::endl;
		}
		for (auto &v: getSortedVertexList()) {
			std::string s1(reinterpret_cast<const char *>(v->getSequence().get()), v->getLength());
			for (auto &edge: outgoingEdgesOf(v)) {
				std::shared_ptr<V> target = getEdgeTarget(edge);
				std::string s2(reinterpret_cast<const char *>(target->getSequence().get()), target->getLength());
				outfile1 << "    " << s1 << " -> " << s2 << "[label=\"" << edge->getMultiplicity() << "\"";
				if (edge->getIsRef())
					outfile1 << ",color=Red";
				outfile1 << "];" << std::endl;
			}
		}
		outfile1 << "}" << std::endl;
		outfile1.close();
	}

	int inDegreeOf(const std::shared_ptr<V> &vertex) {
		return getEdgeContainer(vertex).incoming.size();
	}

	int outDegreeOf(const std::shared_ptr<V> &vector) {
		return getEdgeContainer(vector).outgoing.size();
	}

	phmap::flat_hash_set<std::shared_ptr<E>> &incomingEdgesOf(const std::shared_ptr<V> &vertex) {
		return getEdgeContainer(vertex).getUnmodifiableIncomingEdges();
	}


	phmap::flat_hash_set<std::shared_ptr<E>> &outgoingEdgesOf(const std::shared_ptr<V> &vertex) {
		return getEdgeContainer(vertex).getUnmodifiableOutgoingEdges();
	}

	void removeEdgeFromTouchingVertices(const std::shared_ptr<E> &e) {
		std::shared_ptr<V> source = getEdgeSource(e);
		std::shared_ptr<V> target = getEdgeTarget(e);

		getEdgeContainer(source).removeOutgoingEdge(e);
		getEdgeContainer(target).removeIncomingEdge(e);
	}

	bool addEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex,
	             const std::shared_ptr<E> &e) {
		if (containsEdge(e) || !containsVertex(sourceVertex) || !containsVertex(targetVertex))
			return false;

		if (!allowingMultipleEdges && edgeMap.find(getEdge(sourceVertex, targetVertex)) != edgeMap.end())
			return false;

		if (!allowingLoops && sourceVertex == targetVertex) {
			throw std::invalid_argument(LOOPS_NOT_ALLOWED);
		}

		edgeMap.insert(std::make_pair(e, IntrusiveEdge<V>(sourceVertex, targetVertex)));
		EdgeSet.insert(e);
		getEdgeContainer(sourceVertex).addOutgoingEdge(e);
		getEdgeContainer(targetVertex).addIncomingEdge(e);
		return true;
	}

	bool assertVertexExist(const std::shared_ptr<V> &v) {
		if (vertexMapDirected.find(v) != vertexMapDirected.end())
			return true;
		throw std::invalid_argument("no such vertex in graph.");
	}

	bool containsVertex(const std::shared_ptr<V> &v) {
		return vertexMapDirected.find(v) != vertexMapDirected.end();
	}

	bool containsEdge(const std::shared_ptr<E> &e) {
		return edgeMap.find(e) != edgeMap.end();
	}

	bool isSource(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		return inDegreeOf(v) == 0;
	}

	bool isSink(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		return outDegreeOf(v) == 0;
	}

	std::shared_ptr<uint8_t[]> getAdditionalSequence(const std::shared_ptr<V> &v, bool isSource) {
		return v->getAdditionalSequence(isSource);
	}

	int getAdditionalSequenceLength(const std::shared_ptr<V> &v, bool isSource) {
		return v->getAdditionalSequenceLength(isSource);
	}

	bool removeEdge(const std::shared_ptr<E> &e) {
		if (containsEdge(e)) {
			removeEdgeFromTouchingVertices(e);
			edgeMap.erase(e);
			EdgeSet.erase(e);
			return true;
		}
		return false;
	}

	std::shared_ptr<E> removeEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		std::shared_ptr<E> e = getEdge(sourceVertex, targetVertex);
		if (e != nullptr) {
			removeEdgeFromTouchingVertices(e);
			EdgeSet.erase(e);
			edgeMap.erase(e);
		}
		return e;
	}

	std::shared_ptr<E> addEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		assertVertexExist(sourceVertex);
		assertVertexExist(targetVertex);
		if (!allowingMultipleEdges && containsEdge(sourceVertex, targetVertex))
			return nullptr;

		if (!allowingLoops && sourceVertex == targetVertex)
			throw std::invalid_argument(LOOPS_NOT_ALLOWED);

		std::shared_ptr<E> e = createEdge(sourceVertex, targetVertex);
		if (containsEdge(e))
			return nullptr;
		edgeMap.insert(std::make_pair(e, IntrusiveEdge<V>(sourceVertex, targetVertex)));
		EdgeSet.insert(e);
		getEdgeContainer(sourceVertex).addOutgoingEdge(e);
		getEdgeContainer(targetVertex).addIncomingEdge(e);
		return e;
	}

	bool containsEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		return getEdge(sourceVertex, targetVertex) != nullptr;
	}

	virtual std::shared_ptr<E>
	createEdge(const std::shared_ptr<V> &sourceVertex, const std::shared_ptr<V> &targetVertex) {
		return std::make_shared<E>();
	}

	virtual bool removeVertex(const std::shared_ptr<V> &v) {
		if (containsVertex(v)) {
			removeAllEdges(edgesof(v));
			vertexMapDirected.erase(v);
			VertexSet.erase(v);
			return true;
		}
		return false;
	}

	bool removeAllEdges(const std::vector<std::shared_ptr<E>> &edges) {
		bool modified = false;
		for (const auto &e: edges) {
			modified |= removeEdge(e);
		}
		return modified;
	}

	bool removeAllEdges(const std::list<std::shared_ptr<E>> &edges) {
		bool modified = false;
		for (const auto &e: edges) {
			modified |= removeEdge(e);
		}
		return modified;
	}

	bool removeAllEdges(phmap::flat_hash_set<std::shared_ptr<E>> edges) {
		bool modified = false;
		for (const auto &e: edges) {
			modified |= removeEdge(e);
		}
		return modified;
	}

	bool removeAllVertices(const std::vector<std::shared_ptr<V>> &vertices) {
		bool modified = false;
		for (const std::shared_ptr<V> &v: vertices) {
			modified |= removeVertex(v);
		}
		return modified;
	}

	bool removeAllVertices(const std::list<std::shared_ptr<V>> &vertices) {
		bool modified = false;
		for (const std::shared_ptr<V> &v: vertices) {
			modified |= removeVertex(v);
		}
		return modified;
	}

	bool removeAllVertices(phmap::flat_hash_set<std::shared_ptr<V>> vertices) {
		bool modified = false;
		for (std::shared_ptr<V> v: vertices) {
			modified |= removeVertex(v);
		}
		return modified;
	}

	bool isRefSink(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");

		for (const std::shared_ptr<E> &e: outgoingEdgesOf(v)) {
			if (e->getIsRef())
				return false;
		}

		for (const std::shared_ptr<E> &e: incomingEdgesOf(v)) {
			if (e->getIsRef())
				return true;
		}

		return VertexSet.size() == 1;
	}

	std::shared_ptr<E> incomingEdgeOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		phmap::flat_hash_set<std::shared_ptr<E>> &edgesSet = incomingEdgesOf(v);
		if (edgesSet.size() > 1) {
			throw std::invalid_argument("Cannot get a single incoming edge for a vertex with multiple incoming edges");
		}
		return edgesSet.empty() ? nullptr : *edgesSet.begin();
	}

	std::shared_ptr<E> outgoingEdgeOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		phmap::flat_hash_set<std::shared_ptr<E>> &edgesSet = outgoingEdgesOf(v);
		if (edgesSet.size() > 1) {
			throw std::invalid_argument("Cannot get a single incoming edge for a vertex with multiple incoming edges");
		}
		return edgesSet.empty() ? nullptr : *edgesSet.begin();
	}

	std::shared_ptr<V>
	getNextReferenceVertex(const std::shared_ptr<V> &v, bool allowNonRefPaths, std::shared_ptr<E> blacklistedEdge) {
		if (v == nullptr)
			return nullptr;

		for (const std::shared_ptr<E> &edgeToTest: outgoingEdgesOf(v)) {
			if (edgeToTest->getIsRef()) {
				return getEdgeTarget(edgeToTest);
			}
		}

		if (!allowNonRefPaths)
			return nullptr;

		std::vector<std::shared_ptr<E>> edges;
		for (const std::shared_ptr<E> &edgeToTest: outgoingEdgesOf(v)) {
			if (edgeToTest != blacklistedEdge) {
				edges.template emplace_back(edgeToTest);
			}
			if (edges.size() > 2)
				break;
		}
		return edges.size() == 1 ? getEdgeTarget(edges.at(0)) : nullptr;
	}

	std::shared_ptr<V> getPrevReferenceVertex(const std::shared_ptr<V> &v) {
		if (v == nullptr)
			return nullptr;
		for (const std::shared_ptr<E> &edge: incomingEdgesOf(v)) {
			if (edge->getIsRef())
				return getEdgeSource(edge);
		}
		return nullptr;
	}

	bool isReferenceNode(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		phmap::flat_hash_set<std::shared_ptr<E>> edges = edgesof(v);
		for (const std::shared_ptr<E> &edge: edges) {
			if (edge->getIsRef())
				return true;
		}
		return VertexSet.size() == 1;
	}

	std::shared_ptr<V> getReferenceSourceVertex() {
		for (const std::shared_ptr<V> &vertex: VertexSet) {
			if (isRefSource(vertex))
				return vertex;
		}
		return nullptr;
	}

	std::shared_ptr<V> getReferenceSinkVertex() {
		for (const std::shared_ptr<V> &vertex: VertexSet) {
			if (isRefSink(vertex))
				return vertex;
		}
		return nullptr;
	}


	/**
	 * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
	 *
	 * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
	 * as it requires vertices to not only be connected by a series of directed edges but also prunes away
	 * paths that do not also meet eventually with the reference sink vertex
	 */
	void removePathsNotConnectedToRef(unsigned long size) {
		if (getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr) {
			throw std::invalid_argument("Graph must have ref source and sink vertices");
		}
		std::vector<std::shared_ptr<V>> onPathFromRefSource;
		BaseGraphIterator<V, E> sourceIter = BaseGraphIterator<V, E>(this, getReferenceSourceVertex(), false, true);
		while (sourceIter.hasNext()) {
			onPathFromRefSource.push_back(sourceIter.next());
		}

		phmap::flat_hash_set<std::shared_ptr<V>> onPathFromRefSink;
		BaseGraphIterator<V, E> sinkIter = BaseGraphIterator<V, E>(this, getReferenceSinkVertex(), true, false);
		onPathFromRefSink.reserve(size);
		while (sinkIter.hasNext()) {
			onPathFromRefSink.insert(sinkIter.next());
		}

		phmap::flat_hash_set<std::shared_ptr<V>> verticesToRemove = getVertexSet();
		for (typename std::vector<std::shared_ptr<V>>::iterator iter = onPathFromRefSource.begin();
		     iter != onPathFromRefSource.end(); iter++) {
			if (onPathFromRefSink.find(*iter) != onPathFromRefSink.end())
				verticesToRemove.erase(*iter);
		}

		removeAllVertices(verticesToRemove);

		if (getSinks().size() > 1)
			throw std::length_error("Should have eliminated all but the reference sink");

		if (getSources().size() > 1)
			throw std::length_error("hould have eliminated all but the reference source");
	}

	phmap::flat_hash_set<std::shared_ptr<V>> incomingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		phmap::flat_hash_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: incomingEdgesOf(v)) {
			ret.insert(getEdgeSource(e));
		}
		return ret;
	}

	phmap::flat_hash_set<std::shared_ptr<V>> outgoingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
	    phmap::flat_hash_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: outgoingEdgesOf(v)) {
			ret.insert(getEdgeTarget(e));
		}
		return ret;
	}

	std::vector<std::shared_ptr<V>> vecIncomingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::vector<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: incomingEdgesOf(v)) {
			ret.push_back(getEdgeSource(e));
		}
		return ret;
	}

	std::vector<std::shared_ptr<V>> vecOutgoingVerticesOf(const std::shared_ptr<V> &v) {
		if (v.get() == nullptr)
			throw std::invalid_argument("Attempting to test a null vertex.");
		std::vector<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<E> &e: outgoingEdgesOf(v)) {
			ret.push_back(getEdgeTarget(e));
		}
		return ret;
	}

	phmap::flat_hash_set<std::shared_ptr<V>> getSinks() {
		phmap::flat_hash_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<V> &v: VertexSet) {
			if (isSink(v))
				ret.insert(v);
		}
		return ret;
	}

	phmap::flat_hash_set<std::shared_ptr<V>> getSources() {
		phmap::flat_hash_set<std::shared_ptr<V>> ret;
		for (const std::shared_ptr<V> &v: VertexSet) {
			if (isSource(v))
				ret.insert(v);
		}
		return ret;
	}

	void cleanNonRefPaths() {
		if (getReferenceSourceVertex() == nullptr || getReferenceSinkVertex() == nullptr)
			return;

		phmap::flat_hash_set<std::shared_ptr<E>> edgesToCheck = incomingEdgesOf(getReferenceSourceVertex());
		std::shared_ptr<E> e;

		while (!edgesToCheck.empty()) {
			e = *(edgesToCheck.begin());
			if (!e->getIsRef()) {
				for (const std::shared_ptr<E> &e: incomingEdgesOf(getEdgeSource(e))) {
					edgesToCheck.insert(e);
				}
				removeEdge(e);
			}
			edgesToCheck.erase(e);
		}

		edgesToCheck = outgoingEdgesOf(getReferenceSinkVertex());
		while (!edgesToCheck.empty()) {
			e = *(edgesToCheck.begin());
			if (!e->getIsRef()) {
				for (const std::shared_ptr<E> &e: outgoingEdgesOf(getEdgeTarget(e))) {
					edgesToCheck.insert(e);
				}
				removeEdge(e);
			}
			edgesToCheck.erase(e);
		}

		Specifics<V, E>::removeSingletonOrphanVertices();
	}

	bool isRefSource(const std::shared_ptr<V> &v) {
		return Specifics<V, E>::isRefSource(v);
	}

	void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
        phmap::flat_hash_set<std::shared_ptr<V>> toRemove = VertexSet;
		std::shared_ptr<V> refV = getReferenceSourceVertex();
		if (refV != nullptr) {
			BaseGraphIterator<V, E> iter = BaseGraphIterator<V, E>(this, refV, true, true);
			while (iter.hasNext()) {
				toRemove.erase(iter.next());
			}
		}
		removeAllVertices(toRemove);
	}

	void
	addOrUpdateEdge(const std::shared_ptr<V> &source, const std::shared_ptr<V> &target, const std::shared_ptr<E> &e) {
		if (source.get() == nullptr)
			throw std::invalid_argument("source");
		if (target.get() == nullptr)
			throw std::invalid_argument("target");
		if (e.get() == nullptr)
			throw std::invalid_argument("e");

		std::shared_ptr<E> prev = getEdge(source, target);
		if (prev != nullptr) {
			prev->add(*e);
		} else {
			addEdge(source, target, e);
		}
	}

	phmap::flat_hash_set<std::shared_ptr<E>> getEdgeSet() {
		return EdgeSet;
	}

	bool containsAllVertices(phmap::flat_hash_set<std::shared_ptr<V>> &vertices) {
		if (vertices.empty())
			throw std::invalid_argument("null vertex");
		for (std::shared_ptr<V> v: vertices) {
			if (!containsVertex(v))
				return false;
		}
		return true;
	}

	virtual void reserveSpace(int size) {
		edgeMap.reserve(2 * size);
		EdgeSet.reserve(2 * size);
		vertexMapDirected.reserve(size);
		VertexSet.reserve(size);
	}

};


#endif //MUTECT2CPP_MASTER_DIRECTEDSPECIFICS_H
