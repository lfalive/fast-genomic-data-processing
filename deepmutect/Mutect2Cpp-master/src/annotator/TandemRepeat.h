//
// Created by lhh on 6/15/22.
//

#ifndef MUTECT2CPP_MASTER_TANDEMREPEAT_H
#define MUTECT2CPP_MASTER_TANDEMREPEAT_H

#include "InfoFieldAnnotation.h"

class TandemRepeat : public InfoFieldAnnotation{
public:
    std::vector<std::string> getKeyNames();

    std::shared_ptr<std::map<std::string, AttributeValue>> annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

    static std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>> getNumTandemRepeatUnits(const std::shared_ptr<VariantContext> & vc, const std::shared_ptr<uint8_t[]> & refBasesStartingAtVCWithPad, int refLen, int& len);

    static std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>> getNumTandemRepeatUnits(uint8_t * refBases, int refLen, uint8_t * altBases, int altLen, uint8_t * remainingRefContext, int remainingLen, int& repeatUnitLength);

    static int findRepeatedSubstring(uint8_t * bases, int basesLen);

    static int findNumberOfRepetitions(uint8_t *unit, int length, uint8_t *bases, int len, bool b);
};


#endif //MUTECT2CPP_MASTER_TANDEMREPEAT_H
