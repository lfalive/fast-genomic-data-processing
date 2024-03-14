/**
 * An individual piece of recalibration data. Each bin counts up the number of observations and the number
 * of reference mismatches seen for that combination of covariates.
 */

#ifndef RECAL_DATUM_H
#define RECAL_DATUM_H

#include <limits>


class RecalDatum{
private:
    constexpr static double MULTIPLIER = 100000.0;  //---why is constexpr needed?
    constexpr static double UNINITIALIZED = -1.0;
    constexpr static int SMOOTHING_CONSTANT = 1;
    constexpr static double RESOLUTION_BINS_PER_QUAL = 1.0;


    /**
     * Quals above this value should be capped down to this value (because they are too high)
     * in the base quality score recalibrator
     */
    constexpr static int MAX_GATK_USABLE_Q_SCORE = 40;

    constexpr static long MAX_NUMBER_OF_OBSERVATIONS = std::numeric_limits<int>::max() - 1;

    /**
    * number of bases seen in total
    */
    long numObservations;

    /**
     * Number of bases seen that didn't match the reference
     * (actually sum of the error weights - so not necessarily a whole number)
     * Stored with an internal multiplier to keep it closer to the floating-point sweet spot and avoid numerical error
     * (see https://github.com/broadinstitute/gatk/wiki/Numerical-errors ).
     * However, the value of the multiplier influences the results.
     * For example, you get different results for 1000.0 and 10000.0
     * See MathUtilsUnitTest.testAddDoubles for a demonstration.
     * The value of the MULTIPLIER that we found to give consistent results insensitive to sorting is 10000.0;
    */     //---modified by lhh, remove MULTIPLIER for performance
    double numMismatches;

    /**
     * estimated reported quality score based on combined data's individual q-reporteds and number of observations
     */
    double estimatedQReported;

    /**
     * the empirical quality for datums that have been collapsed together (by read group and reported quality, for example)
     */
    double empiricalQuality;



    /**
     * Calculate and cache the empirical quality score from mismatches and observations (expensive operation)
     */
    void calcEmpiricalQuality(double conditionalPrior);

    static double log10QempPrior(double Qempirical, double Qreported);

    static double log10QempLikelihood(double Qempirical, long nObservations, long nErrors);

public:
    constexpr static char MAX_RECALIBRATION_Q_SCORE = 93;

    static double * log10QempPriorCache;

    RecalDatum(long _numObservations, double _numMismatches, char reportedQuality);

    RecalDatum(RecalDatum & copy);

    // added by lhh, used to calculate log10QempPriorCache array
    static double * CalcLog10QempCache();

    long getNumObservation();

    double getNumMismatch();

    double getEstimatedQReported();

    // calculate the empirical quality every time it is needed
    double getEmpiricalQuality();

    double getEmpiricalQuality(double conditionalPrior);

    // calculate EmpiricalQuality based on estimatedQReported
    void calEmpiricalQuality();

    double GetEmpiricalQuality();

    void increment( double incMismatches);

    /**
     * Add in all of the data from other into this object, updating the reported quality from the expected
     * error rate implied by the two reported qualities
     *
     * @param other  RecalDatum to combine
     */
    void combine(RecalDatum & other);

    /**
     * calculate the expected number of errors given the estimated Q reported and the number of observations
     * in this datum.
     *
     * @return a positive (potentially fractional) estimate of the number of errors
     */
    double calcExpectedErrors();

    void setNumMismatches(double numMismatches);

    void setEmpiricalQuality(double empiricalQuality);

    void setEstimatedQReported(double estimatedQReported);

    static double bayesianEstimateOfEmpiricalQuality(long nObservations, long nErrors, double QReported);


};



#endif