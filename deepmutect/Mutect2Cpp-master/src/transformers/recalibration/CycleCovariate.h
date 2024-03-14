/**
 *  Created by lhh, 2021.3.16
 * */

#ifndef CYCLE_COVARIATE_H
#define CYCLE_COVARIATE_H

#include <string>
#include "htslib/sam.h"
#include "ContextCovariate.h"

class CycleCovariate: public Covariate{
private:
    const static int CUSHION_FOR_INDELS = 4;

    int MaximumCycleValue;
public:
    CycleCovariate(BaseArgument & baseArgument);

    CycleCovariate(RecalibrationArgumentCollection & RAC);

    void recordValues(bam1_t *read, sam_hdr_t *readsHeader, Key * keys);

    /**
      * Computes the encoded value of CycleCovariate's key for the given position at the read.
      * Uses keyFromCycle to do the encoding.
      * @param baseNumber index of the base to compute the key for
      * @param read the read
      *                 (this method throws UserException if the computed absolute value of the cycle number is higher than this value).
      */
    //static int cycleKey(int baseNumber, bam1_t * read);   //---optimized out

    /**
     * Encodes the cycle number as a key.
     */
    int keyFromCycle(int cycle);

    int keyFromValue(string value);

    string getName();

    int maximumKeyValue();

    static string formatKey(int key);
};


#endif
