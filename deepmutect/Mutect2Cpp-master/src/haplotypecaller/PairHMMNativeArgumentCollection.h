//
// Created by lhh on 5/6/22.
//

#ifndef MUTECT2CPP_MASTER_PAIRHMMNATIVEARGUMENTCOLLECTION_H
#define MUTECT2CPP_MASTER_PAIRHMMNATIVEARGUMENTCOLLECTION_H

class PairHMMNativeArgumentCollection{
public:
    int pairHmmNativeThreads = 4;
    bool useDoublePrecision = false;

    PairHMMNativeArgumentCollection()
    {
        pairHmmNativeThreads = 4;
        useDoublePrecision = false;
    }
};

#endif //MUTECT2CPP_MASTER_PAIRHMMNATIVEARGUMENTCOLLECTION_H
