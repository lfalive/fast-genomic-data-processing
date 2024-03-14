//
// Created by 梦想家xixi on 2021/11/12.
//

#ifndef MUTECT2CPP_MASTER_SWPARAMETERS_H
#define MUTECT2CPP_MASTER_SWPARAMETERS_H


class SWParameters {
private:
    int matchValue;
    int mismatchPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;

public:
    SWParameters(int matchValue, int mismatchPenalty, int gapOpenPenalty, int gapExtendPenalty);
    int getGapExtendPenalty() const {return gapExtendPenalty;}
    int getMatchValue() const {return matchValue;}
    int getMismatchPenalty() const {return mismatchPenalty;}
    int getGapOpenPenalty() const {return gapOpenPenalty;}
};


#endif //MUTECT2CPP_MASTER_SWPARAMETERS_H
