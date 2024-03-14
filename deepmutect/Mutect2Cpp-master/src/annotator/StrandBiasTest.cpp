//
// Created by lhh on 6/7/22.
//

#include "StrandBiasTest.h"

vector<vector<int>>
StrandBiasTest::getContingencyTable(AlleleLikelihoods<SAMRecord, Allele> *likelihoods, shared_ptr<VariantContext> vc,
                                    int minCount, vector<string> samples) {
    vector<vector<int>> result;
    if( likelihoods == nullptr || vc == nullptr) {
        return result;
    }

    auto ref = vc->getReference();
    auto allAlts = vc->getAlternateAlleles();

    vector<vector<int>> table(ARRAY_DIM, vector<int>(ARRAY_DIM, 0));
    for(string sample : samples)
    {
        int * sampleTable = new int[ARRAY_SIZE];
        for(int i = 0; i < ARRAY_SIZE; i++) {
            sampleTable[i] = 0;
        }
        auto bestAlleles = likelihoods->bestAllelesBreakingTies(sample);
        for(auto ba : *bestAlleles)
        {
            if(ba->isInformative())
                updateTable(sampleTable, ba->allele, ba->evidence, ref, allAlts);
        }

        if(passesMinimumThreshold(sampleTable, minCount)){
            copyToMainTable(sampleTable, table);
        }

        delete[] sampleTable;
    }
    return table;
}

void
StrandBiasTest::updateTable(int *table, shared_ptr<Allele> allele, shared_ptr<SAMRecord> read, shared_ptr<Allele> ref,
                            std::vector<std::shared_ptr<Allele>> &allAlts) {
    bool matchesRef = allele->equals(*ref, true);
    bool matchesAnyAlt = false;
    for(auto a : allAlts)
    {
        if(*a == *allele)
            matchesAnyAlt = true;
    }

    if( matchesRef || matchesAnyAlt){
        int offset = matchesRef ? 0 : ARRAY_DIM;

        // a normal read with an actual strand
        bool isForward = !read->isReverseStrand();
        table[offset + (isForward ? 0 : 1)]++;
    }
}

bool StrandBiasTest::passesMinimumThreshold(int *data, int minCount) {
    // the ref and alt totals must be greater than MIN_COUNT
    return data[0] + data[1] + data[2] + data[3] > minCount;
}

void StrandBiasTest::copyToMainTable(int *perSampleTable, vector<vector<int>> &mainTable) {
    mainTable[0][0] += perSampleTable[0];
    mainTable[0][1] += perSampleTable[1];
    mainTable[1][0] += perSampleTable[2];
    mainTable[1][1] += perSampleTable[3];
}