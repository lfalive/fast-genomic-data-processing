//
// Created by 梦想家xixi on 2021/11/23.
//

#include "KBestHaplotypeFinder.h"
#include "BaseGraph/DFS_CycleDetect.h"
#include "parallel_hashmap/phmap.h"
#include <memory>
#include <queue>
#include <utility>

KBestHaplotypeFinder::KBestHaplotypeFinder(const std::shared_ptr<SeqGraph> &graph,
                                           phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sources,
                                           phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sinks) : graph(graph) {
	Mutect2Utils::validateArg(graph.get(), "graph cannot be null");
	Mutect2Utils::validateArg(!sources.empty(), "sources cannot be null");
	Mutect2Utils::validateArg(!sinks.empty(), "sinks cannot be null");
	Mutect2Utils::validateArg(graph->containsAllVertices(sources), "source does not belong to the graph");
	Mutect2Utils::validateArg(graph->containsAllVertices(sinks), "sink does not belong to the graph");

	this->graph = DFS_CycleDetect<SeqVertex, BaseEdge>(graph.get()).detectCycles()
	              ? removeCyclesAndVerticesThatDontLeadToSinks(graph, sources, sinks) : graph;
	this->sources = sources;
	this->sinks = sinks;
}

std::shared_ptr<SeqGraph>
KBestHaplotypeFinder::removeCyclesAndVerticesThatDontLeadToSinks(const std::shared_ptr<SeqGraph> &original,
                                                                 phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sources,
                                                                 phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sinks) {
	phmap::flat_hash_set<std::shared_ptr<BaseEdge>> edgesToRemove;
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> vertexToRemove;

	bool foundSomePath = false;
	for (const auto &source: sources) {
		phmap::flat_hash_set<std::shared_ptr<SeqVertex>> parentVertices;
		foundSomePath = findGuiltyVerticesAndEdgesToRemoveCycles(original, source, sinks, edgesToRemove, vertexToRemove,
		                                                         parentVertices) || foundSomePath;
	}
	Mutect2Utils::validateArg(foundSomePath,
	                          "could not find any path from the source vertex to the sink vertex after removing cycles");
	Mutect2Utils::validateArg(!(edgesToRemove.empty() && vertexToRemove.empty()),
	                          "cannot find a way to remove the cycles");

	std::shared_ptr<SeqGraph> result = std::shared_ptr<SeqGraph>(original->clone());

	result->removeAllEdges(edgesToRemove);
	result->removeAllVertices(vertexToRemove);
	return result;
}

bool KBestHaplotypeFinder::findGuiltyVerticesAndEdgesToRemoveCycles(const std::shared_ptr<SeqGraph> &graph,
                                                                    const std::shared_ptr<SeqVertex> &currentVertex,
                                                                    phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sinks,
                                                                    phmap::flat_hash_set<std::shared_ptr<BaseEdge>> &edgesToRemove,
                                                                    phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &verticesToRemove,
                                                                    phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &parentVertices) {
	if (sinks.find(currentVertex) != sinks.end()) {
		return true;
	}
	phmap::flat_hash_set<std::shared_ptr<BaseEdge>> outgoingEdges = graph->outgoingEdgesOf(currentVertex);
	parentVertices.insert(currentVertex);

	bool reachesSink = false;
	for (const auto &edge: outgoingEdges) {
		std::shared_ptr<SeqVertex> child = graph->getEdgeTarget(edge);
		if (parentVertices.find(child) != parentVertices.end()) {
			edgesToRemove.insert(edge);
		} else {
			reachesSink = reachesSink ||
			              findGuiltyVerticesAndEdgesToRemoveCycles(graph, child, sinks, edgesToRemove, verticesToRemove,
			                                                       parentVertices);
		}
	}
	if (!reachesSink) {
		verticesToRemove.insert(currentVertex);
	}
	return reachesSink;
}

KBestHaplotypeFinder::KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, const std::shared_ptr<SeqVertex> &source,
                                           const std::shared_ptr<SeqVertex> &sink) : graph(std::move(graph)) {
	sources.insert(source);
	sinks.insert(sink);
}

KBestHaplotypeFinder::KBestHaplotypeFinder(const std::shared_ptr<SeqGraph> &graph) : graph(graph),
                                                                                     sources(graph->getSources()),
                                                                                     sinks(graph->getSinks()) {}

std::vector<std::shared_ptr<KBestHaplotype>> KBestHaplotypeFinder::findBestHaplotypes(int maxNumberOfHaplotypes) {
	std::vector<std::shared_ptr<KBestHaplotype>> result;
	std::priority_queue<std::shared_ptr<KBestHaplotype>, std::vector<std::shared_ptr<KBestHaplotype>>, KBestHaplotypeComp> queue;
	for (const auto &source: sources) {
		queue.push(std::make_shared<KBestHaplotype>(source, graph));
	}
	phmap::flat_hash_map<std::shared_ptr<SeqVertex>, int> vertexCounts;
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> vertexSet = graph->getVertexSet();
	vertexCounts.reserve(vertexSet.size());
	for (const auto &v: vertexSet) {
		vertexCounts.insert(std::make_pair(v, 0));
	}
	while (!queue.empty() && result.size() < maxNumberOfHaplotypes) {
		std::shared_ptr<KBestHaplotype> pathToExtend = queue.top();
		queue.pop();
		/*int len;
		std::string tmp((char *) pathToExtend->getBases(len).get());
		std::cout.precision(10);
		std::cout.flags(std::ostream::fixed);
		std::cout << pathToExtend->getScore() << "\t" << len << std::endl;
		std::cout << tmp.substr(0, len) << std::endl;*/
		std::shared_ptr<SeqVertex> vertexToExtend = pathToExtend->getLastVertex();
		if (sinks.find(vertexToExtend) != sinks.end()) {
			result.emplace_back(pathToExtend);
		} else {
			if (vertexCounts[vertexToExtend]++ < maxNumberOfHaplotypes) {
				phmap::flat_hash_set<std::shared_ptr<BaseEdge>> outgoingEdges = graph->outgoingEdgesOf(vertexToExtend);
				int totalOutgoingMultiplicity = 0;
				for (const auto &edge: outgoingEdges) {
					totalOutgoingMultiplicity += edge->getMultiplicity();
				}
				for (const auto &edge: outgoingEdges) {
					queue.push(std::make_shared<KBestHaplotype>(pathToExtend, edge, totalOutgoingMultiplicity));
				}
			}
		}
	}
	return result;
}
