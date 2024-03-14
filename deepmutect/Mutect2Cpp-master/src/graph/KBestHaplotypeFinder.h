//
// Created by 梦想家xixi on 2021/11/23.
//

#ifndef MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
#define MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H


#include "SeqGraph.h"
#include "KBestHaplotype.h"
#include <vector>

struct KBestHaplotypeComp {
	bool operator()(const std::shared_ptr<KBestHaplotype> &a, const std::shared_ptr<KBestHaplotype> &b) {
		double aScore = a->getScore(), bScore = b->getScore();
		if (aScore - bScore > 0.0000000001)
			return false;
		if (bScore - aScore > 0.0000000001)
			return true;
		int len1, len2;
		std::shared_ptr<uint8_t[]> base1 = a->getBases(len1), base2 = b->getBases(len2);
		if (len1 != len2)
			return len1 < len2;
		for (int i = 0; i < len1; ++i) {
			if (base1[i] == base2[i]) continue;
			return base1[i] < base2[i];
		}
		return false;
	}
};

class KBestHaplotypeFinder {
private:
	std::shared_ptr<SeqGraph> graph;
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> sinks;
	phmap::flat_hash_set<std::shared_ptr<SeqVertex>> sources;

	static std::shared_ptr<SeqGraph>
	removeCyclesAndVerticesThatDontLeadToSinks(const std::shared_ptr<SeqGraph> &original,
	                                           phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sources,
	                                           phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sinks);

	static bool findGuiltyVerticesAndEdgesToRemoveCycles(const std::shared_ptr<SeqGraph> &graph,
	                                                     const std::shared_ptr<SeqVertex> &currentVertex,
	                                                     phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sinks,
	                                                     phmap::flat_hash_set<std::shared_ptr<BaseEdge>> &edgesToRemove,
	                                                     phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &verticesToRemove,
	                                                     phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &parentVertices);

public:
	KBestHaplotypeFinder(const std::shared_ptr<SeqGraph> &graph,
	                     phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sources,
	                     phmap::flat_hash_set<std::shared_ptr<SeqVertex>> &sinks);

	KBestHaplotypeFinder(std::shared_ptr<SeqGraph> graph, const std::shared_ptr<SeqVertex> &source,
	                     const std::shared_ptr<SeqVertex> &sink);

	KBestHaplotypeFinder(const std::shared_ptr<SeqGraph> &graph);

	std::vector<std::shared_ptr<KBestHaplotype>> findBestHaplotypes(int maxNumberOfHaplotypes);
};


#endif //MUTECT2CPP_MASTER_KBESTHAPLOTYPEFINDER_H
