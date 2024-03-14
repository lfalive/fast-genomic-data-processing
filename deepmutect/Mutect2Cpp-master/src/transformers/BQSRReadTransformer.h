//
// Created by lhh on 3/31/22.
//

#ifndef MUTECT2CPP_MASTER_BQSRREADTRANSFORMER_H
#define MUTECT2CPP_MASTER_BQSRREADTRANSFORMER_H

#include "recalibration/ContextCovariate.h"
#include "recalibration/CycleCovariate.h"
#include "recalibration/StandardCovariateList.h"
#include "recalibration/RecalibrationTables.h"
#include "recalibration/RecalibrationReport.h"
#include "recalibration/ApplyBQSRArgumentCollection.h"

class BQSRReadTransformer {
private:
    shared_ptr<QuantizationInfo> quantizationInfo;
    shared_ptr<RecalibrationTable> recalibrationTables;
    shared_ptr<StandardCovariateList> covariates;


    int preserveQLessThan;
    double globalQScorePrior;
    bool emitOriginalQuals;

    //These fields are created to avoid redoing these calculations for every read
    int totalCovariateCount;
    int specialCovariateCount;

    vector<RecalDatum*> empiricalQualCovsArgs;
    bool useOriginalBaseQualities;

    int cachedReadLength;
    struct Key * cachedKeys;

public:

    BQSRReadTransformer(char * recal_table, ApplyBQSRArgumentCollection & args);
    BQSRReadTransformer(RecalibrationReport report, ApplyBQSRArgumentCollection & args);
    BQSRReadTransformer(shared_ptr<RecalibrationTable> recalibrationTables, shared_ptr<QuantizationInfo> quantizationInfo, shared_ptr<StandardCovariateList> covariates, ApplyBQSRArgumentCollection & args);
    ~BQSRReadTransformer();


    /**
     * Recalibrates the base qualities of a read
     * <p>
     * It updates the base qualities of the read with the new recalibrated qualities (for all event types)
     * <p>
     * Implements a serial recalibration of the reads using the combinational table.
     * First, we perform a positional recalibration, and then a subsequent dinuc correction.
     * <p>
     * Given the full recalibration table, we perform the following preprocessing steps:
     * <p>
     * - calculate the global quality score shift across all data [DeltaQ]
     * - calculate for each of cycle and dinuc the shift of the quality scores relative to the global shift
     * -- i.e., DeltaQ(dinuc) = Sum(pos) Sum(Qual) Qempirical(pos, qual, dinuc) - Qreported(pos, qual, dinuc) / Npos * Nqual
     * - The final shift equation is:
     * <p>
     * Qrecal = Qreported + DeltaQ + DeltaQ(pos) + DeltaQ(dinuc) + DeltaQ( ... any other covariate ... )
     *
     * @param originalRead the read to recalibrate
     */
    void apply(bam1_t * const originalRead);

    // recalibrated quality is bound between 1 and MAX_QUAL
    char getRecalibrateQual(double recalibratedQualDouble);

    static double hierarchicalBayesianQualityEstimate(double epsilon, RecalDatum* empiricalQualRG, RecalDatum* empiricalQualQS, vector<RecalDatum*>& empiricalQualCovs);
};


#endif //MUTECT2CPP_MASTER_BQSRREADTRANSFORMER_H
