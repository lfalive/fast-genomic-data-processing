//
// Created by 梦想家xixi on 2021/10/20.
//

#ifndef MUTECT2CPP_MASTER_MULTISAMPLEEDGE_H
#define MUTECT2CPP_MASTER_MULTISAMPLEEDGE_H

#include "BaseEdge.h"
#include <queue>

class MultiSampleEdge : public BaseEdge {
private:
    int currentSingleSampleMultiplicity{};
    int singleSampleCapacity;
    std::priority_queue<int> singleSampleMultiplicities;

public:
    MultiSampleEdge(bool isRef, int multiplicity, int singleSampleCapacity);

    MultiSampleEdge() : BaseEdge(BaseEdge()), singleSampleCapacity(0) {}

    void flushSingleSampleMultiplicity();

    void incMultiplicity(int incr) override;

    int getPruningMultiplicity() const;

    int getCurrentSingleSampleMultiplicity() const {return currentSingleSampleMultiplicity;}

    bool operator==(const MultiSampleEdge & other) const;

    bool operator<(const MultiSampleEdge & other) const;

};


#endif //MUTECT2CPP_MASTER_MULTISAMPLEEDGE_H
