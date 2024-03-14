//
// Created by lhh on 5/7/22.
//

#ifndef MUTECT2CPP_MASTER_FRAGMENTUTILS_H
#define MUTECT2CPP_MASTER_FRAGMENTUTILS_H

#include "samtools/SAMRecord.h"
#include "QualityUtils.h"

class FragmentUtils {
public:
    constexpr static double DEFAULT_PCR_SNV_ERROR_RATE = 1e-4;
    const static int DEFAULT_PCR_SNV_ERROR_QUAL;
    const static int HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL;

    const static int MISSING_VALUE = 0;

    /**
     * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
     *
     *  Looks at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
     * are not treated independently.  Sets the qualities of firstRead and secondRead to mimic a merged read or
     * nothing if the algorithm cannot create a meaningful one
     * @param pair two overlapping paired reads
     * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
     * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR. May not be negative.
     * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR. May not be negative.
     */
    static void adjustQualsOfOverlappingPairedFragments(std::pair<std::shared_ptr<SAMRecord>, std::shared_ptr<SAMRecord>>& pair, bool setConflictingToZero, int halfOfPcrSnvQual = HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL, int halfOfPcrIndelQual = MISSING_VALUE);
};


#endif //MUTECT2CPP_MASTER_FRAGMENTUTILS_H
