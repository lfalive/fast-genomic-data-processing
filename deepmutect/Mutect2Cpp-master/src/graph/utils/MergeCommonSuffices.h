//
// Created by 梦想家xixi on 2021/11/20.
//

#ifndef MUTECT2CPP_MASTER_MERGECOMMONSUFFICES_H
#define MUTECT2CPP_MASTER_MERGECOMMONSUFFICES_H

# include "VertexBasedTransformer.h"

class MergeCommonSuffices : public VertexBasedTransformer{
public:
    explicit MergeCommonSuffices(SeqGraph* graph) : VertexBasedTransformer(graph) {}

protected:
    bool tryToTransform(std::shared_ptr<SeqVertex> bottom) override;
};


#endif //MUTECT2CPP_MASTER_MERGECOMMONSUFFICES_H
