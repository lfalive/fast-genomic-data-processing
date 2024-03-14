//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_VERTEXBASEDTRANSFORMER_H
#define MUTECT2CPP_MASTER_VERTEXBASEDTRANSFORMER_H


#include "SeqGraph.h"

class VertexBasedTransformer {
private:
	bool dontModifyGraphEvenIfPossible = false;
	SeqGraph *graph;

public:
	bool getDontModifyGraphEvenIfPossible() const { return dontModifyGraphEvenIfPossible; }

	void setDontModifyGraphEvenIfPossible() { dontModifyGraphEvenIfPossible = true; }

	VertexBasedTransformer(SeqGraph *graph);

	SeqGraph *getGraph() { return graph; }

	/**
	 * Merge until the graph has no vertices that are candidates for merging
	 */
	bool transformUntilComplete();

	virtual bool tryToTransform(std::shared_ptr<SeqVertex> v) = 0;
};


#endif //MUTECT2CPP_MASTER_VERTEXBASEDTRANSFORMER_H
