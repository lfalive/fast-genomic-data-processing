//
// Created by 梦想家xixi on 2021/12/1.
//

#include "ReadErrorCorrector.h"
#include "Mutect2Utils.h"
#include "QualityUtils.h"

ReadErrorCorrector::ReadErrorCorrector(int kmerLength, int maxMismatchesToCorrect,
                                       int maxObservationsForKmerToBeCorrectable, uint8_t qualityOfCorrectedBases,
                                       int minObservationsForKmerToBeSolid, bool trimLowQualityBases,
                                       uint8_t minTailQuality, uint8_t *fullReferenceWithPadding, int refLength) : kmerLength(kmerLength),
                                       maxMismatchesToCorrect(maxMismatchesToCorrect), maxObservationsForKmerToBeCorrectable(maxObservationsForKmerToBeCorrectable), qualityOfCorrectedBases(qualityOfCorrectedBases),
                                       minObservationsForKmerToBeSolid(minObservationsForKmerToBeSolid), trimLowQualityBases(trimLowQualityBases), minTailQuality(minTailQuality){
    Mutect2Utils::validateArg(kmerLength > 0, "kmerLength must be > 0");
    Mutect2Utils::validateArg(maxMismatchesToCorrect > 0, "maxMismatchesToCorrect must be >= 1");
    Mutect2Utils::validateArg(qualityOfCorrectedBases >= 2 && qualityOfCorrectedBases <= QualityUtils::MAX_REASONABLE_Q_SCORE, "qualityOfCorrectedBases must be >= 2 and <= MAX_REASONABLE_Q_SCORE");
    maxHomopolymerLengthInRegion = computeMaxHLen(fullReferenceWithPadding, refLength);
}

int ReadErrorCorrector::computeMaxHLen(uint8_t *fullReferenceWithPadding, const int refLength) {
    Mutect2Utils::validateArg(fullReferenceWithPadding, "null is not allowed there.");
    int leftRun = 1;
    int maxRun = 1;
    for ( int i = 1; i < refLength; i++) {
        if ( fullReferenceWithPadding[i] == fullReferenceWithPadding[i-1] ) {
            leftRun++;
        } else {
            leftRun = 1;
        }
    }
    if (leftRun > maxRun) {
        maxRun = leftRun;
    }
    return maxRun;
}

ReadErrorCorrector::ReadErrorCorrector(int kmerLength, uint8_t minTailQuality, int minObservationsForKmerToBeSolid,
                                       uint8_t *fullReferenceWithPadding, int refLength) : kmerLength(kmerLength),
maxMismatchesToCorrect(MAX_MISMATCHES_TO_CORRECT), maxObservationsForKmerToBeCorrectable(MAX_OBSERVATIONS_FOR_KMER_TO_BE_CORRECTABLE), qualityOfCorrectedBases(QUALITY_OF_CORRECTED_BASES),
minObservationsForKmerToBeSolid(minObservationsForKmerToBeSolid), trimLowQualityBases(TRIM_LOW_QUAL_TAILS), minTailQuality(minTailQuality){
    Mutect2Utils::validateArg(kmerLength > 0, "kmerLength must be > 0");
    Mutect2Utils::validateArg(maxMismatchesToCorrect > 0, "maxMismatchesToCorrect must be >= 1");
    Mutect2Utils::validateArg(qualityOfCorrectedBases >= 2 && qualityOfCorrectedBases <= QualityUtils::MAX_REASONABLE_Q_SCORE, "qualityOfCorrectedBases must be >= 2 and <= MAX_REASONABLE_Q_SCORE");
    maxHomopolymerLengthInRegion = computeMaxHLen(fullReferenceWithPadding, refLength);
}

void ReadErrorCorrector::addReadKmers(const std::shared_ptr<SAMRecord>& read) {
    Mutect2Utils::validateArg(read != nullptr, "null is not allowed there");
    if (DONT_CORRECT_IN_LONG_HOMOPOLYMERS && maxHomopolymerLengthInRegion > MAX_HOMOPOLYMER_THRESHOLD) {
        return;
    }
    std::shared_ptr<uint8_t[]> readBases = read->getBases();
    int baseLength = read->getLength();
    for (int offset = 0; offset <= baseLength-kmerLength; offset++ )  {
        countsByKMer->addKmer(new Kmer(readBases,offset,kmerLength),1);
    }
}

void ReadErrorCorrector::addReadsToKmers(const std::vector<std::shared_ptr<SAMRecord>>& reads) {
    for(const auto& read : reads) {
        addReadKmers(read);
    }
}

//SAMRecord *ReadErrorCorrector::correctRead(SAMRecord *inputRead) {
//    Mutect2Utils::validateArg(inputRead, "null is not allowed there");
//    uint8_t * correctedBases = inputRead->getBases();
//    uint8_t * correctedQuals = inputRead->getBaseQualities();
//    int length = inputRead->getLength();
//
//
//}
