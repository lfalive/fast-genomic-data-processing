/**
 * The implementation of QuantizationInfo class
 */

#include <iostream>
#include <fstream>
#include "QuantizationInfo.h"
#include "MathUtils.h"
#include "RecalUtils.h"

QuantizationInfo::QuantizationInfo(RecalibrationTable * tables, int quantizationLevels)
{
    long * qualHistogram = new long[MAX_PHRED_SCORE + 1];
    memset(qualHistogram, 0, (MAX_PHRED_SCORE + 1)* sizeof(long));

    TwoDimensionArray * qualTable = tables->getQualityScoreTable();
    for (TwoDindex i=0; i<tables->getReadGroupNum(); i++)
    {
        for (TwoDindex j=0; j<tables->getQualDimension(); j++)
        {
            RecalDatum * datum = qualTable->operator[](i).operator[](j);
            if (datum == NULL)
                continue;

            int empiricalQual = MathUtils::fastRound(datum->getEmpiricalQuality());
            qualHistogram[empiricalQual] += datum->getNumObservation();
        }
    }
    empiricalQualCounts = qualHistogram;

    quantizeQualityScores(quantizationLevels);
    this->quantizationLevels = quantizationLevels;
}

QuantizationInfo::QuantizationInfo(char *quantizedQuals, long *empiricalQualCounts, int quantizationLevels) :quantizedQuals(quantizedQuals), empiricalQualCounts(empiricalQualCounts), quantizationLevels(quantizationLevels)
{
}

int QuantizationInfo::calculateQuantizationLevels(char *quantizedQuals, int length)
{
    char lastByte = -1;
    int quantizationLevels = 0;
    for(int i=0; i<length; i++)
    {
        if(quantizedQuals[i] != lastByte)
        {
            quantizationLevels++;
            lastByte = quantizedQuals[i];
        }
    }
    return quantizationLevels;
}

void QuantizationInfo::quantizeQualityScores(int nLevels)
{
    QualQuantizer quantizer(empiricalQualCounts, nLevels, QualityUtils::MIN_USABLE_Q_SCORE);
    quantizedQuals = quantizer.getOriginalToQuantizedMap();
}

QuantizationInfo::~QuantizationInfo()
{
    delete[] empiricalQualCounts;   //---delete will case segment fault
    delete[] quantizedQuals;
}

char *QuantizationInfo::getQuantizedQuals()
{
    return quantizedQuals;
}

int QuantizationInfo::getQuantizationLevels()
{
    return quantizationLevels;
}

void QuantizationInfo::noQuantization()
{
    quantizationLevels = QualityUtils::MAX_SAM_QUAL_SCORE;
    for(int i=0; i<quantizationLevels; i++)
    {
        quantizedQuals[i] = i;
    }
}

/*
GATKReportTable* QuantizationInfo::generateReportTable()
{
    GATKReportTable * quantizedTable = new GATKReportTable(RecalUtils::QUANTIZED_REPORT_TABLE_TITLE, "Quality quantization map", 3, Sorting::SORT_BY_COLUMN);
    quantizedTable->addColumn(RecalUtils::QUALITY_SCORE_COLUMN_NAME, "%d");
    quantizedTable->addColumn(RecalUtils::QUANTIZED_COUNT_COLUMN_NAME, "%d");
    quantizedTable->addColumn(RecalUtils::QUANTIZED_VALUE_COLUMN_NAME, "%d");

    for (int qual = 0; qual <= QualityUtils::MAX_SAM_QUAL_SCORE; qual++)
    {
        quantizedTable->set(qual, RecalUtils::QUALITY_SCORE_COLUMN_NAME, qual);
    }
}	*/

void QuantizationInfo::outputToFile(ofstream &output)
{
    /*
     * Table header:
     * #:GATKTable:nColumns:nRows:(DataType for each column):;
     * #:GATKTable:TableName:Description :;
     * key   colA  colB
     * row1  xxxx  xxxxx
    */

    // write the table definition
    output << GATKReportTable::GATKTABLE_HEADER_PREFIX << ":" << QUANTIZEDTABLE_COLUMN_NUM << ":" << QualityUtils::MAX_SAM_QUAL_SCORE + 1;
    for (int i = 0; i < QUANTIZEDTABLE_COLUMN_NUM; ++i) {
        output << GATKReportTable::SEPARATOR << "%d";
    }
    output << GATKReportTable::ENDLINE << endl;
    output << GATKReportTable::GATKTABLE_HEADER_PREFIX << ":" << RecalUtils::QUANTIZED_REPORT_TABLE_TITLE << ":" << "Quality quantization map" << endl;

    // write the column names
    output << RecalUtils::QUALITY_SCORE_COLUMN_NAME << "\t" << RecalUtils::QUANTIZED_COUNT_COLUMN_NAME << "\t" << RecalUtils::QUANTIZED_VALUE_COLUMN_NAME << endl;

    // write the table body
    for(int qual=0; qual <= QualityUtils::MAX_SAM_QUAL_SCORE; qual++)
    {
        output << qual << "\t\t" << empiricalQualCounts[qual] << "\t\t" <<  (int)quantizedQuals[qual] << endl;  //---the format is not tidy
    }
    output << endl;
}

std::ostream & operator<<(std::ostream &os, const QuantizationInfo& quantizationInfo)
{
    for(int qual=0; qual <= QualityUtils::MAX_SAM_QUAL_SCORE; qual++)
    {
        os << qual << "\t" << quantizationInfo.empiricalQualCounts[qual] << "\t" <<  (int)quantizationInfo.quantizedQuals[qual] << endl;
    }
    os << endl;
}