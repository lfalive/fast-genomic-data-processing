/**
 *  A miscellaneous collection of utilities for working with reads, headers, etc.
 *  Static methods only, please.
 */
#ifndef RECALUTILS_H
#define RECALUTILS_H
#include <string>
#include "htslib/sam.h"
#include "RecalibrationTables.h"

class RecalUtils{
public:
    constexpr static int CLIPPING_GOAL_NOT_REACHED = -1;
    const static int EMPIRICAL_Q_REPORTED_DECIMAL_PLACES = 4;
    const static int NUMBER_ERRORS_DECIMAL_PLACES = 2;
    const static string ARGUMENT_REPORT_TABLE_TITLE;
    const static string QUANTIZED_REPORT_TABLE_TITLE;
    const static string READGROUP_REPORT_TABLE_TITLE;
    const static string QUALITY_SCORE_REPORT_TABLE_TITLE;
    const static string ALL_COVARIATES_REPORT_TABLE_TITLE;

    const static string ESTIMATED_Q_REPORTED_COLUMN_NAME;
    const static string QUALITY_SCORE_COLUMN_NAME;
    const static string QUANTIZED_COUNT_COLUMN_NAME;
    const static string QUANTIZED_VALUE_COLUMN_NAME;
    const static string COVARIATE_NAME_COLUMN_NAME;
    const static string COVARIATE_VALUE_COLUMN_NAME;

    const static string NUMBER_OBSERVATIONS_COLUMN_NAME;
    const static string NUMBER_ERRORS_COLUMN_NAME;

    const static int ARGUMENT_COLUMN_INDEX = 0;
    const static int ARGUMENT_VALUE_COLUMN_INDEX = 1;



    /**
    * A marker to tell which end of the read has been clipped
    */
    enum ClippingTail{
        LEFT_TAIL,
        RIGHT_TAIL
    };

    /**
     * Returns the read coordinate corresponding to the requested reference coordinate for a given alignmentStart/Cigar combination.
     *
     * WARNING: if the requested reference coordinate happens to fall inside or just before a deletion (or skipped region) in the read, this function
     * will return the last read base before the deletion (or skipped region). This function returns a
     * Pair(int readCoord, boolean fallsInsideOrJustBeforeDeletionOrSkippedRegion) so you can choose which readCoordinate to use when faced with
     * a deletion (or skipped region).
     *
     * SUGGESTION: Use getReadCoordinateForReferenceCoordinate(GATKSAMRecord, int, ClippingTail) instead to get a
     * pre-processed result according to normal clipping needs. Or you can use this function and tailor the
     * behavior to your needs.
     *
     * @param alignmentStart alignment start of the cigar to the reference
     * @param cigar cigar with which to compute the offset
     * @param refCoord the requested reference coordinate
     * @return the read coordinate corresponding to the requested reference coordinate. (see warning!)
     */
    static int getReadCoordinateForReferenceCoordinate(int alignmentStart, bam1_t * read, int refCoord, ClippingTail tail, bool allowGoalNotReached);

    static pair<int, bool> getReadCoordinateForReferenceCoordinate2(int alignmentStart, bam1_t * read, int refCoord, bool allowGoalNotReached);

    /**
      * Checks if a read starts with an insertion.
      * @return the element if it's a leading insertion or null otherwise
      */
    static uint32_t readStartsWithInsertion(bam1_t * read);

    static void incrementDatumOrPut2Keys(TwoDimensionArray * table, char qual, double isError, int key0, int key1);

    static void incrementDatumOrPut3Keys(ThreeDimensionArray * table, char qual, double isError, int key0, int key1, int key2);

    /**
     * For debugging
     */
    static void printDatum(TwoDimensionArray * table, int length);

};

#endif