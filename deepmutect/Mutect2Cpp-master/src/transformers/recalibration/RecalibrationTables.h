/**
 * Utility class to facilitate base quality score recalibration.
 * author: lhh 2021.3.20
 */

#ifndef RECALIBRATION_TABLE_H
#define RECALIBRATION_TABLE_H

#define MAX_PHRED_SCORE 93

#include <vector>
#include "RecalDatum.h"
#include "StandardCovariateList.h"
#include "boost/multi_array.hpp"

typedef boost::multi_array<RecalDatum *, 2> TwoDimensionArray;
typedef TwoDimensionArray::index TwoDindex;
typedef boost::multi_array<RecalDatum *, 3> ThreeDimensionArray;
typedef ThreeDimensionArray::index ThreeDindex;
typedef boost::multi_array<RecalDatum *, 4> FourDimensionArray;
typedef FourDimensionArray::index FourDindex;

using namespace std;
class RecalibrationTable{
private:
    int qualDimension;
    int numReadGroups;

    int qualityScoreNumber;  // the dimension of output table
    int covTableDimension;  // the dimension of covariate table in the output file

    vector<RecalDatum *> readGroupTable;   //---Is there a better way to represent it?

    //---using boost library
    TwoDimensionArray * qualityScoreTable;
    vector<ThreeDimensionArray*> additionalTable;

    int CovaraiteMaxKey[2];

    // map from covariate name to covariate table
    map<string, ThreeDimensionArray*> covariateToTable;

public:
    RecalibrationTable(int numReadGroups, StandardCovariateList * covariates);

    ~RecalibrationTable();

    TwoDimensionArray * getQualityScoreTable();

    ThreeDimensionArray * getAdditionalTable(int index);

    vector<ThreeDimensionArray*> & getAdditionalTables();

    ThreeDimensionArray * getTableForCovariate(string & covName);

    vector<RecalDatum *>& getReadGroupTable();

    int getReadGroupNum();

    int getQualDimension();

    int getCovariateMaxKey(int index);

    /**
     * merge another table into this one
     * @param table
     */
    void mergeTable(RecalibrationTable * table);

    void increaseQualityScoreNumber();

    int getQualityScoreNumber();

    void increaseCovTableDimension();

    int getCovTableDimension();
};


#endif