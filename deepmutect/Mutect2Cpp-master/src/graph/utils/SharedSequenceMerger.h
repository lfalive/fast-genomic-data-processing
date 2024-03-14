//
// Created by 梦想家xixi on 2021/11/20.
//

#ifndef MUTECT2CPP_MASTER_SHAREDSEQUENCEMERGER_H
#define MUTECT2CPP_MASTER_SHAREDSEQUENCEMERGER_H


#include "SeqGraph.h"

class SharedSequenceMerger {
public:
	static bool canMerge(SeqGraph *graph, const std::shared_ptr<SeqVertex> &v,
	                     phmap::flat_hash_set<std::shared_ptr<SeqVertex>> incomingVertices);

	static bool merge(SeqGraph *graph, const std::shared_ptr<SeqVertex> &v);
};


#endif //MUTECT2CPP_MASTER_SHAREDSEQUENCEMERGER_H
