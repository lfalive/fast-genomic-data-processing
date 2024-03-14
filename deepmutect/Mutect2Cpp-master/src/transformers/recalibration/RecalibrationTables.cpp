/**
 * The implementation of RecalibrationTables class
 */
#include <iostream>
#include "RecalibrationTables.h"
#include <assert.h>

RecalibrationTable::RecalibrationTable(int numReadGroups, StandardCovariateList * covariates)
{
    this->numReadGroups = numReadGroups;
    this->qualDimension = MAX_PHRED_SCORE + 1;

    this->qualityScoreNumber = 0;
    this->covTableDimension = 0;

    this->readGroupTable.reserve(numReadGroups);
    for (int i=0; i<numReadGroups; i++)
    {
        readGroupTable.push_back(NULL);
    }

    this->qualityScoreTable = new TwoDimensionArray(boost::extents[numReadGroups][this->qualDimension]);

    for (int i = 0; i < 2; ++i) {
        CovaraiteMaxKey[i] = covariates->getMaximumKeyValue(i);
    }

    //---TODO: Maybe the data structure can be changed
    this->additionalTable.reserve(2);
    for (int i=0; i<2; i++)
    {
        ThreeDimensionArray * covariateTable = new ThreeDimensionArray(boost::extents[numReadGroups][qualDimension][covariates->getMaximumKeyValue(i)]);
        additionalTable.push_back(covariateTable);
        covariateToTable.insert(pair<string, ThreeDimensionArray*>(covariates->getCovaariateName(i), covariateTable));
    }

}


RecalibrationTable::~RecalibrationTable()
{
    for (int i = 0; i < readGroupTable.size(); ++i) {
        if (readGroupTable[i])
            delete readGroupTable[i];
    }

    for (TwoDindex index0 = 0; index0 < numReadGroups; ++index0) {
        for (TwoDindex index1 = 0; index1 < qualDimension; ++index1) {

                RecalDatum * thisDatum = qualityScoreTable->operator[](index0).operator[](index1);
                if (thisDatum)
                    delete thisDatum;
        }
    }
    delete this->qualityScoreTable;

    for (int i = 0; i < 2; ++i) {
        ThreeDimensionArray * covariateTable = additionalTable[i];

        int CovariateKey = CovaraiteMaxKey[i];
        for (ThreeDindex index0 = 0; index0 < numReadGroups; ++index0) {
            for (ThreeDindex index1 = 0; index1 < qualDimension; ++index1) {
                for (ThreeDindex index2 = 0; index2 < CovariateKey; ++index2) {

                        RecalDatum * thisDatum = covariateTable->operator[](index0).operator[](index1).operator[](index2);
                        if (thisDatum){
                            delete thisDatum;
                	}
		}
            }
        }

        delete covariateTable;
    }

}

TwoDimensionArray * RecalibrationTable::getQualityScoreTable(){
    return qualityScoreTable;
}

ThreeDimensionArray * RecalibrationTable::getAdditionalTable(int index)
{
    return additionalTable[index];
}

vector<ThreeDimensionArray*> & RecalibrationTable::getAdditionalTables()
{
    return additionalTable;
}

ThreeDimensionArray * RecalibrationTable::getTableForCovariate(string &covName)
{
    return covariateToTable[covName];
}

vector<RecalDatum *>& RecalibrationTable::getReadGroupTable()
{
    vector<RecalDatum *>& table = readGroupTable;
    return table;
}

int RecalibrationTable::getReadGroupNum()
{
    return numReadGroups;
}

int RecalibrationTable::getQualDimension()
{
    return qualDimension;
}

int RecalibrationTable::getCovariateMaxKey(int index)
{
    return CovaraiteMaxKey[index];
}

void RecalibrationTable::mergeTable(RecalibrationTable *table)
{
    // merging qualityScoreTable
    TwoDimensionArray * AnotherQualityTable = table->getQualityScoreTable();

    for (TwoDindex index0 = 0; index0 < numReadGroups; ++index0) {
        for (TwoDindex index1 = 0; index1 < qualDimension; ++index1) {
            RecalDatum ** existingDatum = &AnotherQualityTable->operator[](index0).operator[](index1);
            if (*existingDatum)
            {
                RecalDatum ** thisDatum = &qualityScoreTable->operator[](index0).operator[](index1);
                if (*thisDatum)
                    (*thisDatum)->combine(**existingDatum);
                else{
                    //qualityScoreTable->operator[](index0).operator[](index1) = *existingDatum;
		            //AnotherQualityTable->operator[](index0).operator[](index1) = nullptr;
                    *thisDatum = *existingDatum;
                    *existingDatum = nullptr;   // make it convenient to delete it
                }

            }
        }
    }
    
    // merging additionalTable
    for (int i = 0; i < 2; ++i) {
        ThreeDimensionArray * covariateTable = additionalTable[i];
        ThreeDimensionArray * existingTable = table->getAdditionalTable(i);
        //additionalTable[i] = mergeAdditionalTable(covariateTable, existingTable);

        int CovariateKey = CovaraiteMaxKey[i];
        for (ThreeDindex index0 = 0; index0 < numReadGroups; ++index0) {
            for (ThreeDindex index1 = 0; index1 < qualDimension; ++index1) {
                for (ThreeDindex index2 = 0; index2 < CovariateKey; ++index2) {
                    RecalDatum ** existingDatum = &existingTable->operator[](index0).operator[](index1).operator[](index2);
                    if (*existingDatum)
                    {
                        RecalDatum ** thisDatum = &covariateTable->operator[](index0).operator[](index1).operator[](index2);
                        if (*thisDatum)
                            (*thisDatum)->combine(**existingDatum);
                        else
                        {
                            //covariateTable->operator[](index0).operator[](index1).operator[](index2) = *existingDatum;
                            //existingTable->operator[](index0).operator[](index1).operator[](index2) = nullptr;
                            *thisDatum = *existingDatum;
                            *existingDatum = nullptr;
                        }

                    }
                }
            }
        }


    }
    
}

void RecalibrationTable::increaseQualityScoreNumber()
{
    this->qualityScoreNumber ++;
}

int RecalibrationTable::getQualityScoreNumber()
{
    return this->qualityScoreNumber;
}

void RecalibrationTable::increaseCovTableDimension()
{
    this->covTableDimension ++;
}

int RecalibrationTable::getCovTableDimension()
{
    return this->covTableDimension;
}
