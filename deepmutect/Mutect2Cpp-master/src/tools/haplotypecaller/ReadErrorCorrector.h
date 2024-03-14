//
// Created by 梦想家xixi on 2021/12/1.
//

#ifndef MUTECT2CPP_MASTER_READERRORCORRECTOR_H
#define MUTECT2CPP_MASTER_READERRORCORRECTOR_H


#include "KMerCounter.h"
#include "samtools/SAMRecord.h"
#include <cstdint>

class ReadErrorCorrector {
private:
    KMerCounter *countsByKMer;
    std::map<Kmer*, Kmer*> kmerCorrectionMap;
    std::map<Kmer*, std::pair<int*, uint8_t *>> kmerDifferingBases;
    std::map<Kmer*, std::pair<int, int>> kmerDifferingBasesLength;
    int kmerLength;
    bool trimLowQualityBases;
    uint8_t minTailQuality;
    int maxMismatchesToCorrect;
    uint8_t qualityOfCorrectedBases;
    int maxObservationsForKmerToBeCorrectable;
    int maxHomopolymerLengthInRegion;
    int minObservationsForKmerToBeSolid;

    bool doInplaceErrorCorrection = false;
    int MAX_MISMATCHES_TO_CORRECT = 2;
    uint8_t QUALITY_OF_CORRECTED_BASES = 30;
    int MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE = 1;
    bool TRIM_LOW_QUAL_TAILS = false;
    bool DONT_CORRECT_IN_LONG_HOMOPOLYMERS = false;
    int MAX_HOMOPOLYMER_THRESHOLD = 12;
    int numReadsCorrected = 0;
    int numReadsUncorrected = 0;
    int numBasesCorrected = 0;
    int numSolidKmers = 0;
    int numUncorrectableKmers = 0;
    int numCorrectedKmers = 0;

    static int computeMaxHLen(uint8_t* fullReferenceWithPadding, int refLength);

//    SAMRecord* correctRead(SAMRecord* inputRead);

public:
    ReadErrorCorrector(int kmerLength, int maxMismatchesToCorrect, int maxObservationsForKmerToBeCorrectable, uint8_t qualityOfCorrectedBases, int minObservationsForKmerToBeSolid, bool trimLowQualityBases, uint8_t minTailQuality,
                       uint8_t* fullReferenceWithPadding, int refLength);

    ReadErrorCorrector(int kmerLength, uint8_t minTailQuality, int minObservationsForKmerToBeSolid, uint8_t* fullReferenceWithPadding, int refLength);

    void addReadsToKmers(const std::vector<std::shared_ptr<SAMRecord>>& reads);

protected:
    void addReadKmers(const std::shared_ptr<SAMRecord>& read);
};

#endif //MUTECT2CPP_MASTER_READERRORCORRECTOR_H
