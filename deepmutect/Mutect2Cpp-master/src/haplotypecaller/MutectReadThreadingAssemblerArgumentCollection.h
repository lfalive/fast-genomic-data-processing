//
// Created by 梦想家xixi on 2021/12/15.
//

#ifndef MUTECT2CPP_MASTER_MUTECTREADTHREADINGASSEMBLERARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_MUTECTREADTHREADINGASSEMBLERARGUMENTCOLLECTION_H

#include "ReadThreadingAssemblerArgumentCollection.h"

class MutectReadThreadingAssemblerArgumentCollection : public ReadThreadingAssemblerArgumentCollection{
public:
    bool disableAdaptivePruning = false;
    ReadThreadingAssembler* makeReadThreadingAssembler();
};


#endif //MUTECT2CPP_MASTER_MUTECTREADTHREADINGASSEMBLERARGUMENTCOLLECTION_H
