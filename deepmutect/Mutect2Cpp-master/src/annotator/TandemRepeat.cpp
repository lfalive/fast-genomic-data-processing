//
// Created by lhh on 6/15/22.
//

#include "TandemRepeat.h"
#include "VCFConstants.h"
#include "utils/variant/GATKVariantContextUtils.h"

std::vector<std::string> TandemRepeat::getKeyNames() {
    return {VCFConstants::STR_PRESENT_KEY,
            VCFConstants::REPEAT_UNIT_KEY,
            VCFConstants::REPEATS_PER_ALLELE_KEY};
}

std::shared_ptr<std::map<std::string, AttributeValue>>
TandemRepeat::annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc,
                       AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    //todo::filter
    assert(vc != nullptr);
    if ( !vc->isIndel()) {
        return nullptr;
    }
    std::shared_ptr<std::map<std::string, AttributeValue>> map = std::make_shared<std::map<std::string, AttributeValue>>();
    int tmp;
    int repeatLen;
    auto c = ref->getCache(tmp);
    auto result = getNumTandemRepeatUnits(vc, c, tmp, repeatLen);
    if(result.first.empty()) {
        return nullptr;
    }
    std::shared_ptr<uint8_t[]> repeatUnit = result.second;
    std::vector<int> numUnits = result.first;
    std::string s;
    for(int i = 0; i < repeatLen; i++) {
        s += (char)repeatUnit.get()[i];
    }
    map->insert({VCFConstants::STR_PRESENT_KEY, AttributeValue(true)});
    map->insert({VCFConstants::REPEAT_UNIT_KEY, AttributeValue(s)});
    map->insert({VCFConstants::REPEATS_PER_ALLELE_KEY, AttributeValue{numUnits}});
    return map;
}

std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>> TandemRepeat::getNumTandemRepeatUnits(const std::shared_ptr<VariantContext> & vc, const std::shared_ptr<uint8_t[]> & refBasesStartingAtVCWithPad, int refLen, int& len) {
    auto * refBasesStartingAtVCWithoutPad = new uint8_t[refLen-1];
    memcpy(refBasesStartingAtVCWithoutPad, refBasesStartingAtVCWithPad.get()+1, refLen-1);
    std::shared_ptr<Allele> ref = vc->getReference();
    int refAlleleBasesLen = std::max(0, ref->getLength()-1);
    auto * refAlleleBases = new uint8_t[refAlleleBasesLen];
    memcpy(refAlleleBases, ref->getBases().get()+1, refAlleleBasesLen);
//    cout << "refBasesStartingAtVCWithPad " << refLen << " ";
//    for (int i = 0; i < refLen; ++i) cout << refBasesStartingAtVCWithPad.get()[i];
//    cout << endl;
//    cout << "refBasesStartingAtVCWithoutPad " << refLen-1 << " ";
//    for (int i = 0; i < refLen-1; ++i) cout << refBasesStartingAtVCWithoutPad[i];
//    cout << endl;
//    cout << "ref " << ref->getLength() << " ";
//    for (int i = 0; i < ref->getLength(); ++i) cout << ref->getBases().get()[i];
//    cout << endl;
//    cout << "refAlleleBases " << refAlleleBasesLen << " ";
//    for (int i = 0; i < refAlleleBasesLen; ++i) cout << refAlleleBases[i];
//    cout << endl;

    std::shared_ptr<uint8_t[]> repeatUnit;
    std::vector<int> lengths(2, 0);
    bool flag = true;

    for(const auto & allele : vc->getAlternateAlleles()) {
//        cout << "allele " << allele->getBasesLength() << " ";
//        for (int i = 0; i <  allele->getBasesLength(); ++i) cout << allele->getBases()[i];
//        cout << endl;
        int rightLen;
        // make sure that tmp != nullptr, allocate at least 1 uint8_t.
        int allocate_size = max(allele->getBasesLength()-1, 1);
        auto * tmp = new uint8_t[allocate_size];
        memcpy(tmp, allele->getBases().get()+1, allele->getBasesLength()-1);
//        cout << "tmp " << allele->getBasesLength()-1 << " ";
//        for (int i = 0; i <  allele->getBasesLength()-1; ++i) cout << tmp[i];
//        cout << endl;
        std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>> result = getNumTandemRepeatUnits(refAlleleBases, refAlleleBasesLen, tmp, allele->getBasesLength()-1, refBasesStartingAtVCWithoutPad, refLen-1, rightLen);
        delete[] tmp;

        std::vector<int> repetitionCount = result.first;
        if(repetitionCount[0] == 0 || repetitionCount[1] == 0) {
            delete[] refAlleleBases;
            delete[] refBasesStartingAtVCWithoutPad;
            return {};
        }
        if(flag) {
            lengths[0] += repetitionCount[0];
            flag = false;
        }
        lengths[1] += repetitionCount[1];
        repeatUnit = result.second;
        len = rightLen;
    }
    delete[] refAlleleBases;
    delete[] refBasesStartingAtVCWithoutPad;
    return {lengths, repeatUnit};
}

std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>>
TandemRepeat::getNumTandemRepeatUnits(uint8_t * refBases, int refLen, uint8_t * altBases, int altLen, uint8_t * remainingRefContext, int remainingLen, int& repeatUnitLength) {
    uint8_t * longB;
    int longBLen;
    if(altLen > refLen) {
        longB = altBases;
        longBLen = altLen;
    } else {
        longB = refBases;
        longBLen = refLen;
    }

    repeatUnitLength = findRepeatedSubstring(longB, longBLen);
    auto * repeatUnit = new uint8_t[repeatUnitLength];
    memcpy(repeatUnit, longB, repeatUnitLength);

    std::vector<int> repetitionCount = std::vector<int>(2);
    int repetitionsInRef = findNumberOfRepetitions(repeatUnit, repeatUnitLength,refBases, refLen, true);
    auto * tmp = new uint8_t[refLen+remainingLen];
    memcpy(tmp, refBases, refLen);
    memcpy(tmp+refLen, remainingRefContext, remainingLen);
    repetitionCount[0] = findNumberOfRepetitions(repeatUnit, repeatUnitLength, tmp, refLen+remainingLen, true)-repetitionsInRef;
    delete[] tmp;
    tmp = new uint8_t[altLen+remainingLen];
    memcpy(tmp, altBases, altLen);
    memcpy(tmp+altLen, remainingRefContext, remainingLen);
    repetitionCount[1] = findNumberOfRepetitions(repeatUnit, repeatUnitLength, tmp, altLen+remainingLen, true)-repetitionsInRef;
    delete[] tmp;

    return {repetitionCount, std::shared_ptr<uint8_t[]>(repeatUnit)};
}

int TandemRepeat::findRepeatedSubstring(uint8_t *bases, int basesLen) {
    int repLength;
    for (repLength=1; repLength <=basesLen; repLength++) {
        auto * candidateRepeatUnit = new uint8_t[repLength];
        memcpy(candidateRepeatUnit, bases, repLength);
        bool allBasesMatch = true;
        for (int start = repLength; start < basesLen; start += repLength ) {
            // check that remaining of string is exactly equal to repeat unit
            if(basesLen < start + repLength) {
                allBasesMatch = false;
                break;
            }
            auto * basePiece = new uint8_t[repLength];
            memcpy(basePiece, bases+start, repLength);
            if (memcmp(candidateRepeatUnit, basePiece, repLength) != 0) {
                allBasesMatch = false;
                delete[] basePiece;
                break;
            }
            delete[] basePiece;
        }
        delete[] candidateRepeatUnit;
        if (allBasesMatch)
            return repLength;
    }

    return repLength;
}

int TandemRepeat::findNumberOfRepetitions(uint8_t *repeatUnitFull, int repeatUnitLength, uint8_t *testStringFull, int testStringLength, bool leadingRepeats) {
    if(testStringLength == 0) {
        return 0;
    }
    int offsetInRepeatUnitFull = 0;
    int offsetInTestStringFull = 0;
    int lengthDifference = testStringLength - repeatUnitLength;

    if (leadingRepeats) {
        int numRepeats = 0;
        // look forward on the test string
        for (int start = 0; start <= lengthDifference; start += repeatUnitLength) {
            if(memcmp(testStringFull+start + offsetInTestStringFull, repeatUnitFull + offsetInRepeatUnitFull, repeatUnitLength) == 0) {
                numRepeats++;
            } else {
                return numRepeats;
            }
        }
        return numRepeats;
    } else {
        // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
        int numRepeats = 0;
        // look backward on the test string
        for (int start = lengthDifference; start >= 0; start -= repeatUnitLength) {
            if (memcmp(testStringFull + start + offsetInTestStringFull, repeatUnitFull + offsetInRepeatUnitFull, repeatUnitLength) == 0) {
                numRepeats++;
            } else {
                return numRepeats;
            }
        }
        return numRepeats;
    }
}

