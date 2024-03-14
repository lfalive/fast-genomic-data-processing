//
// Created by lhh on 6/7/22.
//

#ifndef MUTECT2CPP_MASTER_STRANDBIASTEST_H
#define MUTECT2CPP_MASTER_STRANDBIASTEST_H

#include "InfoFieldAnnotation.h"

class StrandBiasTest : public InfoFieldAnnotation {
private:
    const static int ARRAY_DIM = 2;
    const static int ARRAY_SIZE = ARRAY_DIM * ARRAY_DIM;

    static void updateTable(int* table, shared_ptr<Allele> allele, shared_ptr<SAMRecord> read, shared_ptr<Allele> ref, std::vector<std::shared_ptr<Allele>>& allAlts);

    static bool passesMinimumThreshold(int* data, int minCount);

    /**
    * Helper method to copy the per-sample table to the main table
    *
    * @param perSampleTable   per-sample table (single dimension)
    * @param mainTable        main table (two dimensions)
    */
    static void copyToMainTable(int* perSampleTable, vector<vector<int>>& mainTable);

public:
    /**
    Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
    *             fw      rc
    *   allele1   #       #
    *   allele2   #       #
    * @return a 2x2 contingency table
    */
    static vector<vector<int>> getContingencyTable( AlleleLikelihoods<SAMRecord, Allele>* likelihoods, shared_ptr<VariantContext> vc, int minCount, vector<string> samples);
};


#endif //MUTECT2CPP_MASTER_STRANDBIASTEST_H
