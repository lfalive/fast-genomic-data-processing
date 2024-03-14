//
// Transformer used to recalibrate the qualities of base
// Created by lhh on 3/31/22.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include "BQSRReadTransformer.h"
#include "recalibration/GATKReportTable.h"
#include "MathUtils.h"


BQSRReadTransformer::BQSRReadTransformer(char * recal_table, ApplyBQSRArgumentCollection & args) : BQSRReadTransformer(RecalibrationReport(recal_table),  args)
{

}

BQSRReadTransformer::BQSRReadTransformer(RecalibrationReport report, ApplyBQSRArgumentCollection & args) : BQSRReadTransformer(report.getRecalibrationTables(), report.getQuantizationInfo(), report.getCovariates(), args)
{

}

BQSRReadTransformer::BQSRReadTransformer(shared_ptr<RecalibrationTable> recalibrationTables, shared_ptr<QuantizationInfo> quantizationInfo, shared_ptr<StandardCovariateList> covariates, ApplyBQSRArgumentCollection & args) :
    quantizationInfo(quantizationInfo), recalibrationTables(recalibrationTables), covariates(covariates), preserveQLessThan(args.PRESERVE_QSCORES_LESS_THAN), globalQScorePrior(args.globalQScorePrior), emitOriginalQuals(args.emitOriginalQuals), useOriginalBaseQualities(args.useOriginalBaseQualities),
    specialCovariateCount(covariates->numberOfSpecialCovariates()), empiricalQualCovsArgs(specialCovariateCount, nullptr), cachedReadLength(0), cachedKeys(
        nullptr){

    if(args.quantizationLevels == 0)
    {
        this->quantizationInfo->noQuantization();
    } else if(args.quantizationLevels > 0 && args.quantizationLevels != quantizationInfo->getQuantizationLevels()) {
        this->quantizationInfo->quantizeQualityScores(args.quantizationLevels);
    }

}

BQSRReadTransformer::~BQSRReadTransformer()
{
    delete[] cachedKeys;
}

void BQSRReadTransformer::apply(bam1_t * const originalRead)
{
    int readLength = originalRead->core.l_qseq;
    if(readLength > cachedReadLength)
    {
        delete[] cachedKeys;
        cachedKeys = new Key[readLength];
    }

    covariates->recordAllValuesInStorage(originalRead, cachedKeys);

    // the rg key is constant over the whole read, the global deltaQ is too
    int rgKey = 0;

    RecalDatum * empiricalQualRG = recalibrationTables->getReadGroupTable()[rgKey];
    if(!empiricalQualRG)
        return;

    uint8_t * quals = bam_get_qual(originalRead);
    double epsilon = globalQScorePrior > 0.0 ? globalQScorePrior: empiricalQualRG->getEstimatedQReported();

    TwoDimensionArray * qualityScoreTable = recalibrationTables->getQualityScoreTable();
    char * quantizedQuals = quantizationInfo->getQuantizedQuals();

    for(int offset = 0; offset < readLength; offset++)
    {
        if(quals[offset] < preserveQLessThan)
            continue;

        for(RecalDatum* & datum : empiricalQualCovsArgs)
        {
            datum = nullptr;
        }

        RecalDatum * empiricalQualQS = (*qualityScoreTable)[rgKey][quals[offset]];

        for(int i=0; i<specialCovariateCount; i++)
        {
            if(i == 0 && cachedKeys[offset].ContextKey >= 0)
                empiricalQualCovsArgs[i] = (*recalibrationTables->getAdditionalTable(i))[rgKey][quals[offset]][cachedKeys[offset].ContextKey];
            else if(i == 1 && cachedKeys[offset].CycleKey >= 0)
                empiricalQualCovsArgs[i] = (*recalibrationTables->getAdditionalTable(i))[rgKey][quals[offset]][cachedKeys[offset].CycleKey];
        }
        double recalibratedQualDouble = hierarchicalBayesianQualityEstimate(epsilon, empiricalQualRG, empiricalQualQS, empiricalQualCovsArgs);

        char recalibratedQualityScore = quantizedQuals[getRecalibrateQual(recalibratedQualDouble)];

        //---print the cached keys and recalibrated quality score
        //cout << (int)quals[offset] << "\t" << cachedKeys[offset].ContextKey << "\t" << cachedKeys[offset].CycleKey << "\t" << (int)recalibratedQualityScore << "\n";

        quals[offset] = (uint8_t)recalibratedQualityScore;


    }

}

char BQSRReadTransformer::getRecalibrateQual(double recalibratedQualDouble)
{
    return QualityUtils::boundQual(MathUtils::fastRound(recalibratedQualDouble), RecalDatum::MAX_RECALIBRATION_Q_SCORE);
}

double BQSRReadTransformer::hierarchicalBayesianQualityEstimate(double epsilon, RecalDatum *empiricalQualRG,
                                                                RecalDatum *empiricalQualQS,
                                                                vector<RecalDatum *> &empiricalQualCovs){
    double globalDeltaQ = empiricalQualRG == nullptr ? 0.0 : empiricalQualRG->getEmpiricalQuality(epsilon) - epsilon;
    double deltaQReported = empiricalQualQS == nullptr ? 0.0 : empiricalQualQS->getEmpiricalQuality(globalDeltaQ + epsilon) - (globalDeltaQ + epsilon);

    double deltaQCovariates = 0.0;
    double conditionalPrior2 = deltaQReported + globalDeltaQ + epsilon;
    for(RecalDatum* empiricalQualCov : empiricalQualCovs)
    {
        if(empiricalQualCov)
        {
            deltaQCovariates += empiricalQualCov->getEmpiricalQuality(conditionalPrior2) - conditionalPrior2;
        }
    }

    return conditionalPrior2 + deltaQCovariates;
}