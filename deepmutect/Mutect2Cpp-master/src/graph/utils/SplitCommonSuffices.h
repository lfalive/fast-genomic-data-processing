//
// Created by 梦想家xixi on 2021/11/19.
//

#ifndef MUTECT2CPP_MASTER_SPLITCOMMONSUFFICES_H
#define MUTECT2CPP_MASTER_SPLITCOMMONSUFFICES_H


#include "VertexBasedTransformer.h"

class SplitCommonSuffices : public VertexBasedTransformer{
private:
    std::set<std::shared_ptr<SeqVertex>> alreadySplit;

public:
    explicit SplitCommonSuffices(SeqGraph* graph) : VertexBasedTransformer(graph) {}

protected:
    bool tryToTransform(std::shared_ptr<SeqVertex> bottom) override;
};


#endif //MUTECT2CPP_MASTER_SPLITCOMMONSUFFICES_H
