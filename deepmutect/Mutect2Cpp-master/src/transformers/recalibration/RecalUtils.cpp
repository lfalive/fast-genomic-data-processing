/**
 * The implementation of class RecalUtils
 */
#include <iostream>
#include <algorithm>
#include "RecalUtils.h"
#include "RecalDatum.h"

const string RecalUtils::ARGUMENT_REPORT_TABLE_TITLE = "Arguments";
const string RecalUtils::QUANTIZED_REPORT_TABLE_TITLE = "Quantized";
const string RecalUtils::READGROUP_REPORT_TABLE_TITLE = "RecalTable0";
const string RecalUtils::QUALITY_SCORE_REPORT_TABLE_TITLE = "RecalTable1";
const string RecalUtils::ALL_COVARIATES_REPORT_TABLE_TITLE = "RecalTable2";

const string RecalUtils::ESTIMATED_Q_REPORTED_COLUMN_NAME = "EstimatedQReported";
const string RecalUtils::QUALITY_SCORE_COLUMN_NAME = "QualityScore";
const string RecalUtils::QUANTIZED_COUNT_COLUMN_NAME = "Count";
const string RecalUtils::QUANTIZED_VALUE_COLUMN_NAME = "QuantizedScore";

const string RecalUtils::COVARIATE_NAME_COLUMN_NAME = "CovariateName";
const string RecalUtils::COVARIATE_VALUE_COLUMN_NAME = "CovariateValue";

const string RecalUtils::NUMBER_OBSERVATIONS_COLUMN_NAME = "Observations";
const string RecalUtils::NUMBER_ERRORS_COLUMN_NAME = "Errors";


int RecalUtils::getReadCoordinateForReferenceCoordinate(int alignmentStart, bam1_t *read, int refCoord, ClippingTail tail, bool allowGoalNotReached)
{

    pair<int, bool> result = getReadCoordinateForReferenceCoordinate2(alignmentStart, read, refCoord, allowGoalNotReached);
    int readCoord = result.first;

    // Corner case one: clipping the right tail and falls on deletion, move to the next
    // read coordinate. It is not a problem for the left tail because the default answer
    // from getReadCoordinateForReferenceCoordinate is to give the previous read coordinate.
    if (result.second && tail == ClippingTail::RIGHT_TAIL)
        readCoord++;

    uint32_t firstElementIsInsertion = readStartsWithInsertion(read);

    if (readCoord == 0 && tail == ClippingTail::LEFT_TAIL && firstElementIsInsertion != 0)
        readCoord = bam_cigar_oplen(firstElementIsInsertion) < read->core.l_qseq - 1 ? bam_cigar_oplen(firstElementIsInsertion) : read->core.l_qseq - 1;

    return readCoord;
}

pair<int, bool> RecalUtils::getReadCoordinateForReferenceCoordinate2(int alignmentStart, bam1_t *read, int refCoord, bool allowGoalNotReached)
{
    int readBases = 0;
    int refBases = 0;
    bool fallsInsideDeletionOrSkippedRegion = false;
    bool endJustBeforeDeletionOrSkippedRegion = false;
    bool fallsInsideOrJustBeforeDeletionOrSkippedRegion = false;
    uint32_t * cigar = bam_get_cigar(read);

    int goal = refCoord - alignmentStart;
    if (goal < 0)
    {
        if (allowGoalNotReached)
            return pair<int, bool>(CLIPPING_GOAL_NOT_REACHED, false);
        else
            throw "Somehow the requested coordinate is not covered by the read. Too many deletions?";
    }

    bool goalReached = refBases == goal;

    for (unsigned i=0; !goalReached && i<read->core.n_cigar ; i++)
    {
        int shift = 0;
        uint32_t cigarElement = cigar[i];
        int Operator = bam_cigar_op(cigarElement);
        if ( (bam_cigar_type(Operator) & 2) || Operator == BAM_CSOFT_CLIP){
            if (refBases + bam_cigar_oplen(cigarElement) < goal){
                shift = bam_cigar_oplen(cigarElement);
            } else{
                shift = goal - refBases;
            }
            refBases += shift;
        }
        goalReached = refBases == goal;
        if (!goalReached && (bam_cigar_type(Operator) & 1)){
            readBases += bam_cigar_oplen(cigarElement);
        }

        if (goalReached){
            // Is this base's reference position within this cigar element? Or did we use it all?
            bool endsWithinCigar = shift < bam_cigar_oplen(cigarElement);

            // If it isn't, we need to check the next one. There should *ALWAYS* be a next one
            // since we checked if the goal coordinate is within the read length, so this is just a sanity check.
            if (!endsWithinCigar && i == read->core.n_cigar - 1){
                if (allowGoalNotReached)
                    return pair<int, bool>(CLIPPING_GOAL_NOT_REACHED, false);
                else
                    throw "Reference coordinate corresponds to a non-existent base in the read. This should never happen";
            }

            uint32_t nextCigarElement = 0;

            // if we end inside the current cigar element, we just have to check if it is a deletion (or skipped region)
            if(endsWithinCigar){
                fallsInsideDeletionOrSkippedRegion = ( Operator == BAM_CDEL || Operator == BAM_CREF_SKIP);
            }
            else{
                i++;
                nextCigarElement = cigar[i];    //---in Java, here is cigarElementIterator.next()

                // if it's an insertion, we need to clip the whole insertion before looking at the next element
                if (bam_cigar_op(nextCigarElement) == BAM_CINS)
                {
                    readBases += bam_cigar_oplen(nextCigarElement);
                    if (i == read->core.n_cigar - 1)
                    {
                        if (allowGoalNotReached)
                            return pair<int, bool>(CLIPPING_GOAL_NOT_REACHED, false);
                        else
                            throw "Reference coordinate corresponds to a non-existent base in the read. This should never happen";
                    }

                    i++;
                    nextCigarElement = cigar[i];
                }

                // if it's a deletion (or skipped region), we will pass the information on to be handled downstream.
                endJustBeforeDeletionOrSkippedRegion = (bam_cigar_op(nextCigarElement) == BAM_CDEL || bam_cigar_op(nextCigarElement) == BAM_CREF_SKIP);
            }

            fallsInsideOrJustBeforeDeletionOrSkippedRegion = endJustBeforeDeletionOrSkippedRegion || fallsInsideDeletionOrSkippedRegion;

            // If we reached our goal outside a deletion (or skipped region), add the shift
            if (!fallsInsideOrJustBeforeDeletionOrSkippedRegion && (bam_cigar_type(Operator) & 1)){
                readBases += shift;
            }
                // If we reached our goal just before a deletion (or skipped region) we need
                // to add the shift of the current cigar element but go back to it's last element to return the last
                // base before the deletion (or skipped region) (see warning in function contracts)
            else if (endJustBeforeDeletionOrSkippedRegion && (bam_cigar_type(Operator) & 1)){
                readBases += shift - 1;
            }// If we reached our goal inside a deletion (or skipped region), or just between a deletion and a skipped region,
                // then we must backtrack to the last base before the deletion (or skipped region)
            else if (fallsInsideDeletionOrSkippedRegion ||
                     (endJustBeforeDeletionOrSkippedRegion && bam_cigar_op(nextCigarElement) == BAM_CREF_SKIP) ||
                     (endJustBeforeDeletionOrSkippedRegion && bam_cigar_op(nextCigarElement) == BAM_CDEL)){
                readBases--;
            }
        }
    }

    if (!goalReached)
    {
        if (allowGoalNotReached)
            return pair<int, bool>(CLIPPING_GOAL_NOT_REACHED, false);
        else
            throw "Somehow the requested coordinate is not covered by the read. Too many deletions?";
    }

    return pair<int, bool>(readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion);
}

uint32_t RecalUtils::readStartsWithInsertion(bam1_t *read)
{
    uint32_t * cigar = bam_get_cigar(read);
    for (unsigned i=0; i<read->core.n_cigar; i++)
    {
        if (bam_cigar_op(cigar[i]) == BAM_CINS){
            return cigar[i];
        }
        else if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP && bam_cigar_op(cigar[i]) != BAM_CSOFT_CLIP){
            break;
        }
    }
    return 0;
}

void RecalUtils::incrementDatumOrPut2Keys(TwoDimensionArray *table, char qual, double isError, int key0, int key1)
{
    TwoDindex index0 = key0;
    TwoDindex index1 = key1;
    //RecalDatum * existingDatum = table->operator[](index0).operator[](index1); //Maybe using reference is better

    if (table->operator[](index0).operator[](index1) == NULL)
    {
        table->operator[](index0).operator[](index1) = new RecalDatum(1, isError, qual);
    } else
    {
        table->operator[](index0).operator[](index1)->increment(isError);
        //table->operator[](index0).operator[](index1)->getEmpiricalQuality();
    }
}

void RecalUtils::incrementDatumOrPut3Keys(ThreeDimensionArray *table, char qual, double isError, int key0, int key1, int key2)
{
    ThreeDindex index0 = key0;
    ThreeDindex index1 = key1;
    ThreeDindex index2 = key2;
    RecalDatum * existingDatum = table->operator[](index0).operator[](index1).operator[](index2);

    if (existingDatum == NULL)
    {
        table->operator[](index0).operator[](index1).operator[](index2) = new RecalDatum(1, isError, qual);
    }
    else{
        table->operator[](index0).operator[](index1).operator[](index2)->increment(isError);
    }
}

void RecalUtils::printDatum(TwoDimensionArray *table, int length)
{
    TwoDindex index0 = 0;
    TwoDindex index1 ;

    for (index1 = 0; index1 < length; index1++)
    {
        if (table->operator[](index0).operator[](index1) != NULL)
        {
            cout << table->operator[](index0).operator[](index1)->getEstimatedQReported() << "\t";
            cout << table->operator[](index0).operator[](index1)->getNumObservation() << "\t";
            cout << table->operator[](index0).operator[](index1)->getNumMismatch() << "\t";
            cout << table->operator[](index0).operator[](index1)->getEmpiricalQuality() << endl;
        }
    }


}
