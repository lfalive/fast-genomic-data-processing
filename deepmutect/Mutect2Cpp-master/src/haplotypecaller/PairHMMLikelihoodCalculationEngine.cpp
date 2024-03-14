//
// Created by 梦想家xixi on 2022/3/1.
//

#include <cassert>
#include <memory>
#include <algorithm>
#include <utility>
#include "PairHMMLikelihoodCalculationEngine.h"
#include "utils/pairhmm/VectorLoglessPairHMM.h"
#include "MathUtils.h"
#include "ReadUtils.h"
#include "utils/variant/GATKVariantContextUtils.h"
#include "utils/Utils.h"


using namespace std;

double PairHMMLikelihoodCalculationEngine::INITIAL_QSCORE = 40.0;
double PairHMMLikelihoodCalculationEngine::EXPECTED_ERROR_RATE_PER_BASE = 0.02;
char PairHMMLikelihoodCalculationEngine::constantGCP = 10;

template<> double AlleleLikelihoods<SAMRecord, Haplotype>::NATURAL_LOG_INFORMATIVE_THRESHOLD = MathUtils::log10ToLog(LOG_10_INFORMATIVE_THRESHOLD);

PairHMMLikelihoodCalculationEngine::PairHMMLikelihoodCalculationEngine(char constantGCP, PairHMMNativeArgumentCollection& args,
                                                                       double log10globalReadMismappingRate,
                                                                       PCRErrorModel pcrErrorModel,
                                                                       char baseQualityScoreThreshold)
{
    assert(constantGCP >= 0);
    assert(log10globalReadMismappingRate <= 0);
    PairHMMLikelihoodCalculationEngine::constantGCP = constantGCP;
    this->log10globalReadMismappingRate = log10globalReadMismappingRate;
    this->pcrErrorModel = pcrErrorModel;
    this->pairHMM = std::make_shared<VectorLoglessPairHMM>(args);

    initializePCRErrorModel();

    assert(baseQualityScoreThreshold >= QualityUtils::MIN_USABLE_Q_SCORE);
    this->baseQualityScoreThreshold = baseQualityScoreThreshold;
}

PairHMMLikelihoodCalculationEngine::~PairHMMLikelihoodCalculationEngine() = default;

char PairHMMLikelihoodCalculationEngine::getErrorModelAdjustedQual(int repeatLength, double rateFactor)
{
    return (char)std::max(MIN_ADJUSTED_QSCORE, MathUtils::fastRound(INITIAL_QSCORE - exp( repeatLength / (rateFactor * M_PI) )) + 1 );
}

void PairHMMLikelihoodCalculationEngine::initializePCRErrorModel()
{
    if(!hasRateFactor())
        return;

    pcrIndelErrorModelCache = shared_ptr<char[]>(new char[MAX_REPEAT_LENGTH + 1]);
    double rateFactor = getRateFactor();
    for(int i=0; i<=MAX_REPEAT_LENGTH; i++)
    {
        pcrIndelErrorModelCache[i] = getErrorModelAdjustedQual(i, rateFactor);
    }
}

void PairHMMLikelihoodCalculationEngine::computeReadLikelihoods(SampleMatrix<SAMRecord, Haplotype> *likelihoods)
{
    shared_ptr<vector<shared_ptr<SAMRecord>>> processedReads = modifyReadQualities(likelihoods->evidence());
    auto gapContinuationPenalties = buildGapContinuationPenalties(*processedReads, constantGCP);

    // Run the PairHMM to calculate the log10 likelihood of each (processed) reads' arising from each haplotype
    if(BOOST_UNLIKELY(pairHMM->is_use_trietree_optimize))
        pairHMM->computeLog10Likelihoods_trie_unique(likelihoods, *processedReads, gapContinuationPenalties);
    else
        pairHMM->computeLog10Likelihoods(likelihoods, *processedReads, gapContinuationPenalties);

    delete gapContinuationPenalties;
}

AlleleLikelihoods<SAMRecord, Haplotype>* PairHMMLikelihoodCalculationEngine::computeReadLikelihoods(AssemblyResultSet &assemblyResultSet, std::vector<std::string>& samples,
                                                                std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> &perSampleReadList)
{
    std::shared_ptr<std::vector<std::shared_ptr<Haplotype>>> haplotypeList = assemblyResultSet.getHaplotypeList();

    initializePairHMM(*haplotypeList, perSampleReadList);

    // Add likelihoods for each sample's reads to our result
    auto* result = new AlleleLikelihoods<SAMRecord, Haplotype>(samples, haplotypeList, perSampleReadList);

    int sampleCount = samples.size();
    for(int i=0; i<sampleCount; i++)
    {
        computeReadLikelihoods(result->sampleMatrix(i));
    }

    result->normalizeLikelihoods(log10globalReadMismappingRate);
    result->filterPoorlyModeledEvidence(&PairHMMLikelihoodCalculationEngine::log10MinTrueLikelihood, EXPECTED_ERROR_RATE_PER_BASE);

    return result;
}

void PairHMMLikelihoodCalculationEngine::initializePairHMM(vector<shared_ptr<Haplotype>> &haplotypes,
                                                           map<string, vector<shared_ptr<SAMRecord>>> &perSampleReadList)
{
    int readMaxLength = 0;
    for(auto& pair : perSampleReadList)
    {
        for(shared_ptr<SAMRecord> & read : pair.second)
        {
            if(read->getLength() > readMaxLength)
                readMaxLength = read->getLength();
        }
    }

    int haplotypeMaxLength = 0;
    for(auto & haplotype : haplotypes)
    {
        if(haplotype->getLength() > haplotypeMaxLength)
            haplotypeMaxLength = haplotype->getLength();
    }

    // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
    pairHMM->initialize(haplotypes, perSampleReadList, readMaxLength, haplotypeMaxLength);
}

shared_ptr<vector<shared_ptr<SAMRecord>>> PairHMMLikelihoodCalculationEngine::modifyReadQualities(vector<shared_ptr<SAMRecord>>& reads)
{
    auto result = make_shared<vector<shared_ptr<SAMRecord>>>();
    result->reserve(reads.size());
    for(auto & read : reads)
    {
        //---use getBases() or getBasesNoCopy()?
        shared_ptr<uint8_t[]> readBases = read->getBases();

        // NOTE -- must clone anything that gets modified here so we don't screw up future uses of the read
        //Using close here is justified - it's an array of primitives.
        shared_ptr<uint8_t[]> readQuals = read->getBaseQualities();
        int length = 0;
        shared_ptr<uint8_t[]> readInsQuals = ReadUtils::getBaseInsertionQualities(read, length);
        shared_ptr<uint8_t[]> readDelQuals = ReadUtils::getBaseDeletionQualities(read, length);

        applyPCRErrorModel(length, readBases.get(), readInsQuals.get(), readDelQuals.get());
        capMinimumReadQualities(*read, read->getBaseQualitiesLength(), readQuals.get(), readInsQuals.get(), readDelQuals.get(), baseQualityScoreThreshold);

        // Create a new copy of the read and sets its base qualities to the modified versions.
        result->emplace_back(createQualityModifiedRead(*read, read->getLength(), readBases, readQuals, readInsQuals, readDelQuals));
    }

    return result;
}

void PairHMMLikelihoodCalculationEngine::applyPCRErrorModel(int length, uint8_t *readBases, uint8_t *readInsQuals,
                                                            uint8_t *readDelQuals) {
    for(int i=1; i<length; i++)
    {
        int repeatLength = findTandemRepeatUnits(readBases, length, i-1);
        readInsQuals[i-1] = (uint8_t)std::min(0xff & readInsQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        readDelQuals[i-1] = (uint8_t)std::min(0xff & readDelQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
    }
}

bool ArraysEqual(vector<uint8_t> & array1, vector<uint8_t> & array2)
{
    int length1 = array1.size();
    int length2 = array2.size();
    if(length1 != length2)
        return false;

    for(int i=0; i<length1; i++)
    {
        if(array1[i] != array2[i])
            return false;
    }
    return true;
}


int PairHMMLikelihoodCalculationEngine::findTandemRepeatUnits(uint8_t *readBases, int length, int offset)
{
    int maxBW = 0;
    vector<uint8_t> bestBWRepeatUnit = {readBases[offset]};

    for(int str = 1; str <= MAX_STR_UNIT_LENGTH; str++)
    {
        // fix repeat unit length
        //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
        if(offset+1-str < 0)
        {
            break;
        }

        // get backward repeat unit and # repeats
        maxBW = GATKVariantContextUtils::findNumberOfRepetitions(readBases, length, offset - str + 1,  str , readBases, length, 0, offset + 1, false);
        if(maxBW > 1)
        {
            bestBWRepeatUnit.clear();
            bestBWRepeatUnit.reserve(str);
            for(int i=offset-str+1; i<offset+1; i++)
            {
                bestBWRepeatUnit.push_back(readBases[i]);
            }
            break;
        }
    }
    vector<uint8_t> bestRepeatUnit = bestBWRepeatUnit;
    int maxRL = maxBW;

    if(offset < length - 1)
    {
        vector<uint8_t> bestFWRepeatUnit = {readBases[offset+1]};
        int maxFW = 0;

        for(int str = 1; str <= MAX_STR_UNIT_LENGTH; str++)
        {
            // fix repeat unit length
            //edge case: if candidate tandem repeat unit falls beyond edge of read, skip
            if (offset+str+1 > length) {
                break;
            }

            // get forward repeat unit and # repeats
            maxFW = GATKVariantContextUtils::findNumberOfRepetitions(readBases, length, offset + 1, str, readBases, length, offset + 1, length-offset -1, true);
            if (maxFW > 1) {
                //bestFWRepeatUnit = Arrays.copyOfRange(readBases, offset + 1, offset+str+1);
                bestFWRepeatUnit.clear();
                bestFWRepeatUnit.reserve(str);
                for(int i=offset+1; i<offset+str+1; i++)
                {
                    bestFWRepeatUnit.push_back(readBases[i]);
                }
                break;
            }
        }

        // if FW repeat unit = BW repeat unit it means we're in the middle of a tandem repeat - add FW and BW components
        if(ArraysEqual(bestFWRepeatUnit, bestBWRepeatUnit)){
            maxRL = maxFW + maxBW;
            bestRepeatUnit = bestFWRepeatUnit;
        } else {
            // tandem repeat starting forward from current offset.
            // It could be the case that best BW unit was different from FW unit, but that BW still contains FW unit.
            // For example, TTCTT(C) CCC - at (C) place, best BW unit is (TTC)2, best FW unit is (C)3.
            // but correct representation at that place might be (C)4.
            // Hence, if the FW and BW units don't match, check if BW unit can still be a part of FW unit and add
            // representations to total

            maxBW = GATKVariantContextUtils::findNumberOfRepetitions(bestFWRepeatUnit.data(), bestFWRepeatUnit.size(), readBases, offset+1,
                                                                     false);
            maxRL = maxFW + maxBW;
        }
    }

    if(maxRL > MAX_REPEAT_LENGTH) {
        maxRL = MAX_REPEAT_LENGTH;
    }
    return maxRL;
}

void PairHMMLikelihoodCalculationEngine::capMinimumReadQualities(SAMRecord &read, int readQualsLength,
                                                                 uint8_t *readQuals, uint8_t *readInsQuals,
                                                                 uint8_t *readDelQuals,
                                                                 char baseQualityScoreThreshold){
    for(int i=0; i<readQualsLength; i++)
    {
        readQuals[i] = (uint8_t)min(0xff & readQuals[i], read.getMappingQuality());
        readQuals[i] = setToFixedValueIfTooLow(readQuals[i], baseQualityScoreThreshold, QualityUtils::MIN_USABLE_Q_SCORE);
        readInsQuals[i] = setToFixedValueIfTooLow( readInsQuals[i], QualityUtils::MIN_USABLE_Q_SCORE,       QualityUtils::MIN_USABLE_Q_SCORE );
        readDelQuals[i] = setToFixedValueIfTooLow( readDelQuals[i], QualityUtils::MIN_USABLE_Q_SCORE,       QualityUtils::MIN_USABLE_Q_SCORE );
    }
}

uint8_t PairHMMLikelihoodCalculationEngine::setToFixedValueIfTooLow(uint8_t currentVal, uint8_t minQual,
                                                                    uint8_t fixedQual){
    return currentVal < minQual ? fixedQual : currentVal;
}

shared_ptr<SAMRecord> PairHMMLikelihoodCalculationEngine::createQualityModifiedRead(SAMRecord& read, int length, std::shared_ptr<uint8_t[]> readBases, std::shared_ptr<uint8_t[]> baseQualities, std::shared_ptr<uint8_t[]> baseInsertionQualities, std::shared_ptr<uint8_t[]> baseDeletionQualities)
{
    shared_ptr<SAMRecord> processedRead = ReadUtils::emptyRead(read);
    processedRead->setBases(std::move(readBases), length);
    processedRead->setBaseQualities(std::move(baseQualities), length);
    ReadUtils::setInsertionBaseQualities(processedRead, std::move(baseInsertionQualities), length);
    ReadUtils::setDeletionBaseQualities(processedRead, std::move(baseDeletionQualities), length);
    return processedRead;
}

phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>* PairHMMLikelihoodCalculationEngine::buildGapContinuationPenalties(vector<shared_ptr<SAMRecord>>& reads, char gapPenalty)
{
    auto result = new phmap::flat_hash_map<SAMRecord*, shared_ptr<char[]>>(reads.size());
    for(auto& read: reads)
    {
        result->emplace(pair<SAMRecord*, shared_ptr<char[]>>(read.get(), Utils::dupBytes(gapPenalty, read->getLength())));
    }
    return result;
}

double PairHMMLikelihoodCalculationEngine::log10MinTrueLikelihood(shared_ptr<SAMRecord> read, double maximumErrorPerBase)
{
    double maxErrorsForRead = min(2.0, ceil(read->getLength() * maximumErrorPerBase));
    double log10QualPerBase = -4.0;
    return maxErrorsForRead * log10QualPerBase;
}