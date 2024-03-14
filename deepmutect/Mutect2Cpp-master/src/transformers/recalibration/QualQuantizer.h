/**
 * Created by lhh, 2021.3.26
 */

/**
 * A general algorithm for quantizing quality score distributions to use a specific number of levels
 *
 * Takes a histogram of quality scores and a desired number of levels and produces a
 * map from original quality scores -> quantized quality scores.
 */

#ifndef QUAL_QUANTIZER
#define QUAL_QUANTIZER

#include <set>
#include "RecalibrationTables.h"
using namespace std;

class QualQuantizer{

    class QualInterval{
    //protected:
    public:
        int qStart, qEnd, fixedQual, level;
        long nObservations, nErrors;
        set<QualInterval> *  subIntervals;


    //public:     //---protected methods can't be accessed by other class
        int mergeOrder; //---Maybe useless

        int minQual;    //---added by lhh

        QualInterval(int qStart, int qEnd, long nObservation, long nErrors, int level, int fixedQual, int minQual);

        QualInterval(int qStart, int qEnd, long nObservation, long nErrors, int level, int fixedQual, int minQual, set<QualInterval > * subIntervals);

        QualInterval(const QualInterval & interval);

        ~QualInterval();

        /**
         * Compare two QualInterval's qStart
         */
        int compareTo(QualInterval* qualInterval);

        /**
         * Create an interval representing the merge of this interval and toMerge
         * @param toMerge
         * @return newly created merged QualInterval
         */
        QualInterval* merge(QualInterval* toMerge);

        double getPenalty();

        double calcPenalty (double globalErrorRate) const ;

        double getErrorRate() const;

        bool operator< (const QualQuantizer::QualInterval & toCompare) const
        {
            return (this->qStart < toCompare.qStart);
        }

        /**
         * @return the QUAL of the error rate of this interval, or the fixed qual if this interval was created with a fixed qual.
         */
        char getQual() const;

        void freeSubIntervals();
    };


private:

    /**
     * Inputs to the QualQuantizer
     */
    int nLevels;

    long * nObservationsPerQual;

    const static int NQualsInHistogram = MAX_PHRED_SCORE + 1;

    /**
     * Map from original qual(e.g. Q30) to new quantized qual(e.g., Q28).
     * Has the same range as nObservationserQual
     */
     char * originalToQuantizeMap;

     /**
      * Sorted set of qual intervals
      * After quantize() this data structure contains only the top-level qual intervals
      */
      set<QualInterval > * quantizedIntervals;

     int getNQualsInHistogram();

     /**
      * Given a final forest of intervals, constructs a list mapping list.get(i) => quantized qual to use for original quality score i
      *
      * This function should be called only once to initialize the corresponding
      * cached value in this object, as the calculation is a bit costly.
      */
     char * intervalsToMap(set<QualInterval > * intervals);
public:
    int minInterestingQual;

    /**
     * Creates a QualQuantizer for the histogram that has nLevels
     *
     * Note this is the only interface to the system.  After creating this object
     * the map can be obtained via getOriginalToQuantizedMap()
     *
     * @param nObservationsPerQual A histogram of counts of bases with quality scores.  Note that
     *  this histogram must start at 0 (i.e., get(0) => count of Q0 bases) and must include counts all the
     *  way up to the largest quality score possible in the reads.  OK if the histogram includes many 0
     *  count bins, as these are quantized for free.
     * @param nLevels the desired number of distinct quality scores to represent the full original range.  Must
     *  be at least 1.
     * @param minInterestingQual All quality scores <= this value are considered uninteresting and are freely
     *  merged together.  For example, if this value is 10, then Q0-Q10 are all considered free to merge, and
     *  quantized into a single value. For ILMN data with lots of Q2 bases this results in a Q2 bin containing
     *  all data with Q0-Q10.
     */
    QualQuantizer(long * nObservationsPerQual, int nLevels, int minInterestingQual);

    ~QualQuantizer();

    /**
     * Main method for computing the quantization intervals.
     *
     * Invoked in the constructor after all input variables are initialized.  Walks
     * over the inputs and builds the min. penalty forest of intervals with exactly nLevel
     * root nodes.  Finds this min. penalty forest via greedy search, so is not guarenteed
     * to find the optimal combination.
     *
     * TODO: develop a smarter algorithm
     *
     * @return the forest of intervals with size == nLevels
     */
     set<QualQuantizer::QualInterval > * quantize();   //---remember to add QualQuantizer:: because it's nested class

    /**
     * Helper function that finds and merge together the lowest penalty pair of intervals
     */
     void mergeLowestPenaltyIntervals(set<QualInterval> * intervals);

     char * getOriginalToQuantizedMap();
};



#endif