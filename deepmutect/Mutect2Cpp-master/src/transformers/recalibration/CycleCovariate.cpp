/**
 * The implementation of class CycleCovariate
 * */

#include <iostream>
#include <math.h>
#include <string.h>
//#include "SamRead.h"
#include "CycleCovariate.h"

CycleCovariate::CycleCovariate(BaseArgument & baseArgument):MaximumCycleValue(baseArgument.MaximumCycleValue) {

}

CycleCovariate::CycleCovariate(RecalibrationArgumentCollection &RAC) : MaximumCycleValue(RAC.MAXIMUM_CYCLE_VALUE)
{

}

void CycleCovariate::recordValues(bam1_t *read, sam_hdr_t *readsHeader, Key * keys)
{
    int readLength = read->core.l_qseq;
    bool isNegStrand = bam_is_rev(read);
    bool isSecondInPair = isPaired(read) && getSecondOfPairFlagUnchecked(read);

    int readOrderFactor = isSecondInPair ? -1 : 1;
    int increment;
    int cycle;
    if (isNegStrand)
    {
        cycle = readLength * readOrderFactor;
        increment = -1 * readOrderFactor;
    } else{
        cycle = readOrderFactor;
        increment = readOrderFactor;
    }

    for (int i=0; i<readLength; i++)
    {
        int substitutionKey = keyFromCycle(cycle);
        keys[i].CycleKey = substitutionKey;
        cycle += increment;
    }
}
/*
int CycleCovariate::cycleKey(int baseNumber, bam1_t *read) {
    bool isNegStrand = bam_is_rev(read);
    bool isSecondInPair = isPaired(read) && getSecondOfPairFlagUnchecked(read);
    int readLength = read->core.l_qseq;

    int readOrderFactor = isSecondInPair ? -1 : 1;
    int increment;
    int cycle;
    if (isNegStrand)
    {
        cycle = readLength * readOrderFactor;
        increment = -1 * readOrderFactor;
    } else{
        cycle = readOrderFactor;
        increment = readOrderFactor;
    }

    cycle += baseNumber * increment;

    return keyFromCycle(cycle, MAXIMUM_CYCLE_VALUE);
}*/

int CycleCovariate::keyFromCycle(int cycle)
{
    int result = abs(cycle);
    if (result > MaximumCycleValue) {
        throw "a larger cycle than the maximum was detected.  Please use the --maximum_cycle_value argument to increase this value (at the expense of requiring more memory to run)";
    }

    result <<= 1; // shift so we can add the "sign" bit
    if ( cycle < 0 ) {
        result++; // negative cycles get the lower-most bit set
    }
    return result;

}

int CycleCovariate::keyFromValue(string value)
{
    return keyFromCycle(stoi(value));
}

string CycleCovariate::getName()
{
    return "Cycle";
}

int CycleCovariate::maximumKeyValue()
{
    return (MaximumCycleValue << 1) + 1;
}

string CycleCovariate::formatKey(int key)
{
    int cycle = key >> 1; // shift so we can remove the "sign" bit
    if ( (key & 1) != 0 ) { // is the last bit set?
        cycle *= -1; // then the cycle is negative
    }
    return to_string(cycle);
}