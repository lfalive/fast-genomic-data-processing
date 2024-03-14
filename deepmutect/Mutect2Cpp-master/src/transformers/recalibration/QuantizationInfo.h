/**
 * Class that encapsulate the information necessary for quality score quantization for BQSR
 */

#ifndef QUANTIZATION_INFO_H
#define QUANTIZATION_INFO_H

#include <vector>
#include "RecalibrationTables.h"
#include "QualQuantizer.h"
#include "QualityUtils.h"
#include "GATKReportTable.h"

#define QUANTIZEDTABLE_COLUMN_NUM 3

using namespace std;

class QuantizationInfo{
private:
    //vector<char> quantizedQuals;
    long * empiricalQualCounts;
    char * quantizedQuals;
    int quantizationLevels;



public:
    QuantizationInfo(char * quantizedQuals, long * empiricalQualCounts, int quantizationLevels);

    QuantizationInfo(RecalibrationTable * tables, int quantizationLevels);

    ~QuantizationInfo();

    void quantizeQualityScores(int nLevels);

    char * getQuantizedQuals();

    int getQuantizationLevels();

    GATKReportTable* generateReportTable();

    void outputToFile(ofstream & output);

    /**
     * calculate the quantizationLevels
     * @param quantizedQuals
     * @param length        the length of the array
     * @return
     */
    static int calculateQuantizationLevels(char * quantizedQuals, int length);

    void noQuantization();

    friend std::ostream & operator<<(std::ostream &os, const QuantizationInfo& quantizationInfo);
};

#endif