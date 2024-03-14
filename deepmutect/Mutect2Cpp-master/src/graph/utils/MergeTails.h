//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_MERGETAILS_H
#define MUTECT2CPP_MASTER_MERGETAILS_H


#include "VertexBasedTransformer.h"
#include "SharedVertexSequenceSplitter.h"

class MergeTails : public VertexBasedTransformer{
public:
    static const int MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES = 10;

    explicit MergeTails(SeqGraph* graph) : VertexBasedTransformer(graph){}

protected:
    bool tryToTransform(std::shared_ptr<SeqVertex> top) override;
};


#endif //MUTECT2CPP_MASTER_MERGETAILS_H
