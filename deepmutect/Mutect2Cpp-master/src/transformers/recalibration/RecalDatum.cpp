/**
 * The implementation of RecalDatum class
 *
 */

#include <iostream>
#include <limits>
#include <math.h>
#include "assert.h"
#include "RecalDatum.h"
#include "QualityUtils.h"
#include "MathUtils.h"

double * RecalDatum::log10QempPriorCache = RecalDatum::CalcLog10QempCache();    //static initialization should be outside any function

/*
RecalDatum::RecalDatum(long _numObservations, double _numMismatches, char reportedQuality)
{
    if (_numObservations < 0)
        throw "numObservations < 0";
    if (_numMismatches < 0.0)
        throw "numMismatch < 0";
    if (reportedQuality < 0)
        throw "reportedQuality < 0";

    this->numObservations = _numObservations;
    //this->numMismatches = (_numMismatches * MULTIPLIER);
    this->numMismatches = _numMismatches;
    this->estimatedQReported = reportedQuality;
    this->empiricalQuality = UNINITIALIZED;
    //this->empiricalQuality = estimatedQReported;
}
*/

RecalDatum::RecalDatum(long _numObservations, double _numMismatches, char reportedQuality) : numObservations(_numObservations), numMismatches(_numMismatches),
                                    estimatedQReported(reportedQuality), empiricalQuality(UNINITIALIZED)
{
    assert(estimatedQReported >= 0);
    assert(numMismatches >= 0);

    //calEmpiricalQuality();
}

RecalDatum::RecalDatum(RecalDatum &copy)
{
    this->numObservations = copy.numObservations;
    this->numMismatches = copy.numMismatches;
    this->estimatedQReported = copy.estimatedQReported;
    this->empiricalQuality = copy.empiricalQuality;
}

// f(x) = a*exp(-((x - b)^2 / (2*c^2)))
// Note that a is the height of the curve's peak, b is the position of the center of the peak, and c controls the width of the "bell".
double * RecalDatum::CalcLog10QempCache()
{
    double GF_a = 0.9;
    double GF_b = 0.0;
    double GF_c = 0.5;

    double temp = 2 * GF_c * GF_c;
    double * log10QempPriorCache = new double[MAX_GATK_USABLE_Q_SCORE + 1];
    for (int i = 0; i <= MAX_GATK_USABLE_Q_SCORE; ++i) {
        double value = GF_a * exp(-(i - GF_b) * (i - GF_b) / temp); //---to be tested
        double log10Prior = log10(value);
        if ( isinf(log10Prior))
            log10Prior = -std::numeric_limits<double>::max();
        log10QempPriorCache[i] = log10Prior;
    }

    return log10QempPriorCache;
}

long RecalDatum::getNumObservation()
{
    return numObservations;
}

double RecalDatum::getNumMismatch()
{
    //return this->numMismatches/MULTIPLIER;
    return this->numMismatches;
}

double RecalDatum::getEstimatedQReported()
{
    return estimatedQReported;
}

double RecalDatum::GetEmpiricalQuality()
{
    return empiricalQuality;
}

double RecalDatum::getEmpiricalQuality()
{
    return getEmpiricalQuality(getEstimatedQReported());
}

double RecalDatum::getEmpiricalQuality(double conditionalPrior)
{
    if (empiricalQuality == UNINITIALIZED)
    {
        calcEmpiricalQuality(conditionalPrior);
    }
    return empiricalQuality;
}

void RecalDatum::calEmpiricalQuality()
{
    if (empiricalQuality == UNINITIALIZED)
    {
        calcEmpiricalQuality(estimatedQReported);
    }
}

void RecalDatum::calcEmpiricalQuality(double conditionalPrior)
{
    long mismatches = (long)(getNumMismatch() + 0.5) + SMOOTHING_CONSTANT;
    long observations = getNumObservation() + SMOOTHING_CONSTANT + SMOOTHING_CONSTANT;

    double empiricalQual = RecalDatum::bayesianEstimateOfEmpiricalQuality(observations, mismatches, conditionalPrior);

    empiricalQuality = empiricalQual < MAX_RECALIBRATION_Q_SCORE ? empiricalQual : (double)MAX_RECALIBRATION_Q_SCORE;
}

double RecalDatum::bayesianEstimateOfEmpiricalQuality(long nObservations, long nErrors, double QReported)
{
    int numBins = (QualityUtils::MAX_REASONABLE_Q_SCORE + 1) * (int)RESOLUTION_BINS_PER_QUAL;

    int MLEbin = 0;
    double * log10Posteriors = new double[numBins];
    for (int i=0; i<numBins; i++)
    {
        double QEmpOfBin = i / RESOLUTION_BINS_PER_QUAL;
        log10Posteriors[i] = log10QempPrior(QEmpOfBin, QReported) + log10QempLikelihood(QEmpOfBin, nObservations, nErrors);

        if (log10Posteriors[i] > log10Posteriors[MLEbin])
            MLEbin = i;
    }
    delete[] log10Posteriors;
    return MLEbin / RESOLUTION_BINS_PER_QUAL;

}

double RecalDatum::log10QempPrior(double Qempirical, double Qreported)
{
    int difference = std::min(abs( (int)(Qempirical - Qreported)), MAX_GATK_USABLE_Q_SCORE);
    return log10QempPriorCache[difference];
}

double RecalDatum::log10QempLikelihood(double Qempirical, long nObservations, long nErrors)
{
    if (nObservations == 0)
        return 0.0;

    // the binomial code requires ints as input (because it does caching).  This should theoretically be fine because
    // there is plenty of precision in 2^31 observations, but we need to make sure that we don't have overflow
    // before casting down to an int.
    if (nObservations > MAX_NUMBER_OF_OBSERVATIONS)
    {
        // we need to decrease nErrors by the same fraction that we are decreasing nObservations
        double fraction = (double) MAX_NUMBER_OF_OBSERVATIONS / (double)nObservations;
        nErrors = round((double) nErrors * fraction);
        nObservations = MAX_NUMBER_OF_OBSERVATIONS;
    }

    //this is just a straight binomial PDF
    double log10Prob = MathUtils::log10BinomialProbability((int)nObservations, (int)nErrors, Qempirical / (-10.0));
    if (isinf(log10Prob) || isnan(log10Prob))
        log10Prob = -std::numeric_limits<double>::max();

    return log10Prob;
}

void RecalDatum::increment( double incMismatches)
{
    numObservations++;
    //numObservations += incObservations;
    //numMismatches += (incMismatches * MULTIPLIER);  //---multiplication is too expensive!
    numMismatches += incMismatches;

}

void RecalDatum::combine(RecalDatum &other)
{
    double sumErrors = this->calcExpectedErrors() + other.calcExpectedErrors();
    //increment(other.numObservations, other.getNumMismatch());
    numObservations += other.numObservations;
    numMismatches += other.getNumMismatch();
    estimatedQReported = -10 * log10(sumErrors / numObservations);

}

double RecalDatum::calcExpectedErrors() {
    return numObservations * pow(10.0, estimatedQReported / -10.0);
}

void RecalDatum::setNumMismatches(double numMismatches)
{
    if (numMismatches < 0.0)
        throw "numMismatches < 0";
    //this->numMismatches = (numMismatches * MULTIPLIER);
    this->numMismatches = numMismatches;
    empiricalQuality = UNINITIALIZED;
}

void RecalDatum::setEmpiricalQuality(double empiricalQuality)
{
    if (empiricalQuality < 0.0 || isinf(empiricalQuality) || isnan(empiricalQuality))
        throw "empiricalQuality is not qualified";
    this->empiricalQuality = empiricalQuality;
}

void RecalDatum::setEstimatedQReported(double estimatedQReported)
{
    if (estimatedQReported < 0.0 || isinf(estimatedQReported) || isnan(estimatedQReported))
        throw "empiricalQuality is not qualified";
    this->estimatedQReported = estimatedQReported;
    empiricalQuality = UNINITIALIZED;
    //calcExpectedErrors();
}
