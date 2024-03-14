//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_SEQGRAPH_H
#define MUTECT2CPP_MASTER_SEQGRAPH_H

#include "BaseGraph/DirectedSpecifics.h"
#include "SeqVertex.h"
#include "BaseEdge.h"
#include <list>

class SeqGraph : public DirectedSpecifics<SeqVertex, BaseEdge> {
private:
	int kmerSize;

private:
	/**
	 * How many cycles of the graph simplifications algorithms will we run before
	 * thinking something has gone wrong and throw an exception?
	 */
	static const int MAX_REASONABLE_SIMPLIFICATION_CYCLES = 100;

	/**
	 * Is source vertex potentially a start of a linear chain of vertices?
	 *
	 * We are a start of a zip chain if our out degree is 1 and either the
	 * the vertex has no incoming connections or 2 or more (we must start a chain) or
	 * we have exactly one incoming vertex and that one has out-degree > 1 (i.e., source's incoming
	 * vertex couldn't be a start itself
	 *
	 * @param source a non-null vertex
	 * @return true if source might start a linear chain
	 */
	bool isLinearChainStart(const std::shared_ptr<SeqVertex> &source);

	/**
	 * Get all of the vertices in a linear chain of vertices starting at zipStart
	 *
	 * Build a list of vertices (in order) starting from zipStart such that each sequential pair of vertices
	 * in the chain A and B can be zipped together.
	 *
	 * @param zipStart a vertex that starts a linear chain
	 * @return a list of vertices that comprise a linear chain starting with zipStart.  The resulting
	 *         list will always contain at least zipStart as the first element.
	 */
	std::list<std::shared_ptr<SeqVertex>> traceLinearChain(const std::shared_ptr<SeqVertex> &zipStart);

	bool mergeLinearChain(std::list<std::shared_ptr<SeqVertex>> &linearChain);

	static std::shared_ptr<SeqVertex> mergeLinearChainVertices(std::list<std::shared_ptr<SeqVertex>> &vertices);

	/**
	 * Run one full cycle of the graph simplification algorithms
	 * @return true if any algorithms said they did some simplification
	 */
	bool simplifyGraphOnce(int iteration);

public:
	int getKmerSize() const { return kmerSize; }

	std::shared_ptr<BaseEdge>
	createEdge(std::shared_ptr<SeqVertex> sourceVertex, std::shared_ptr<SeqVertex> targetVertrx);

	SeqGraph(int kmer) : kmerSize(kmer), DirectedSpecifics<SeqVertex, BaseEdge>() {}

	SeqGraph(SeqGraph &seqGraph);

	/**
	 * Zip up all of the simple linear chains present in this graph.
	 *
	 * Merges together all pairs of vertices in the graph v1 -> v2 into a single vertex v' containing v1 + v2 sequence
	 *
	 * Only works on vertices where v1's only outgoing edge is to v2 and v2's only incoming edge is from v1.
	 *
	 * If such a pair of vertices is found, they are merged and the graph is update.  Otherwise nothing is changed.
	 *
	 * @return true if any such pair of vertices could be found, false otherwise
	 */
	bool zipLinearChains();

	/**
	* Simplify this graph, merging vertices together and restructuring the graph in an
	* effort to minimize the number of overall vertices in the graph without changing
	* in any way the sequences implied by a complex enumeration of all paths through the graph.
	*/
	void simplifyGraph();

	void simplifyGraph(int maxCycles);

	SeqGraph *clone();

	void printGraphSize(const std::string &info);
};


#endif //MUTECT2CPP_MASTER_SEQGRAPH_H
