//
// Created by 梦想家xixi on 2021/11/24.
//

#ifndef MUTECT2CPP_MASTER_KBESTHAPLOTYPE_H
#define MUTECT2CPP_MASTER_KBESTHAPLOTYPE_H

#include "path/Path.h"
#include "SeqGraph.h"
#include "haplotype/Haplotype.h"

class KBestHaplotype : public Path<SeqVertex, BaseEdge>{
private:
    double score;
    bool isReference = false;

public:
    double getScore() {return score;}
    bool getIsReference() {return isReference;}
    KBestHaplotype(std::shared_ptr<SeqVertex> initialVertex, std::shared_ptr<DirectedSpecifics<SeqVertex, BaseEdge>> graph);
    KBestHaplotype(const std::shared_ptr<KBestHaplotype>& p, const std::shared_ptr<BaseEdge>& edge, int totalOutgoingMultiplicity);
    std::shared_ptr<Haplotype> getHaplotype();
};


#endif //MUTECT2CPP_MASTER_KBESTHAPLOTYPE_H
