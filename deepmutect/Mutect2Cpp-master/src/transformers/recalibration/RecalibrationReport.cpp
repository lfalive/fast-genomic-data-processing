//
// Created by lhh on 4/1/22.
//

#include <iostream>
#include <memory>
#include <cassert>
#include "RecalibrationReport.h"
#include "RecalUtils.h"

RecalibrationReport::RecalibrationReport(char *recal_table) : RecalibrationReport(GATKReport(recal_table))
{

}

RecalibrationReport::RecalibrationReport(GATKReport report)
{
    argumentTable = report.getTable(RecalUtils::ARGUMENT_REPORT_TABLE_TITLE);
    RAC = shared_ptr<RecalibrationArgumentCollection>(initializeArgumentCollectionTable(*argumentTable));

    GATKReportTable* quantizedTable = report.getTable(RecalUtils::QUANTIZED_REPORT_TABLE_TITLE).get();
    quantizationInfo = shared_ptr<QuantizationInfo>(initializeQuantizationTable(*quantizedTable));

    covariates = make_shared<StandardCovariateList>(*RAC);
    recalibrationTables = make_shared<RecalibrationTable>(report.getReadGroups().size(), covariates.get());

    parseReadGroupTable(*report.getTable(RecalUtils::READGROUP_REPORT_TABLE_TITLE), recalibrationTables->getReadGroupTable());
    parseQualityScoreTable(*report.getTable(RecalUtils::QUALITY_SCORE_REPORT_TABLE_TITLE), recalibrationTables->getQualityScoreTable());
    parseAllCovariatesTable(*report.getTable(RecalUtils::ALL_COVARIATES_REPORT_TABLE_TITLE), *recalibrationTables);

}

RecalibrationReport::~RecalibrationReport()
{
}

RecalibrationArgumentCollection* RecalibrationReport::initializeArgumentCollectionTable(GATKReportTable &table)
{
    RecalibrationArgumentCollection * RAC = new RecalibrationArgumentCollection();
    for(int i=0; i<table.getNumRows(); i++)
    {
        string argument = table.get(i, RecalUtils::ARGUMENT_COLUMN_INDEX);
        string value = table.get(i, RecalUtils::ARGUMENT_VALUE_COLUMN_INDEX);
        if(argument == "mismatches_context_size")
            RAC->MISMATCHES_CONTEXT_SIZE = stoi(value);
        else if(argument == "indels_context_size")
            RAC->INDELS_CONTEXT_SIZE = stoi(value);
        else if(argument == "mismatches_default_quality")
            RAC->MISMATCHES_DEFAULT_QUALITY = stoi(value);
        else if(argument == "insertions_default_quality")
            RAC->INSERTIONS_DEFAULT_QUALITY = stoi(value);
        else if(argument == "deletions_default_quality")
            RAC->DELETIONS_DEFAULT_QUALITY = stoi(value);
        else if(argument == "maximum_cycle_value")
            RAC->MAXIMUM_CYCLE_VALUE = stoi(value);
        else if(argument == "low_quality_tail")
            RAC->LOW_QUAL_TAIL = stoi(value);
        else if(argument == "quantizing_levels")
            RAC->QUANTIZING_LEVELS = stoi(value);
        else if(argument == "binary_tag_name")
            RAC->BINARY_TAG_NAME = value;
    }
    return RAC;
}

QuantizationInfo* RecalibrationReport::initializeQuantizationTable(GATKReportTable &table)
{
    char* quals = new char[QualityUtils::MAX_SAM_QUAL_SCORE + 1];
    long* counts = new long[QualityUtils::MAX_SAM_QUAL_SCORE + 1];
    for(int i=0; i<table.getNumRows(); i++)
    {
        quals[i] = stoi(table.get(i, 2));
        counts[i] = stol(table.get(i, 1));
    }
    return new QuantizationInfo(quals, counts, QuantizationInfo::calculateQuantizationLevels(quals, QualityUtils::MAX_SAM_QUAL_SCORE + 1));
}

void RecalibrationReport::parseReadGroupTable(GATKReportTable& reportTable, vector<RecalDatum *> &readGroupTable)
{
    for(int i=0; i<reportTable.getNumRows(); i++)
    {
        int ReadGroupNum = i;
        readGroupTable[i] = getRecalDatum(reportTable, ReadGroupNum, true);
    }
}

void RecalibrationReport::parseQualityScoreTable(GATKReportTable &reportTable, TwoDimensionArray *qualTable)
{
    int tempQUALarray[2];
    for(int i=0; i<reportTable.getNumRows(); i++)
    {
        tempQUALarray[0] = 0;   // TODO: how to make it more elegant?
        tempQUALarray[1] = stoi(reportTable.get(i, RecalUtils::QUALITY_SCORE_COLUMN_NAME));

        (*qualTable)[tempQUALarray[0]][tempQUALarray[1]] = getRecalDatum(reportTable, i, false);
    }
}

void RecalibrationReport::parseAllCovariatesTable(GATKReportTable &reportTable, RecalibrationTable &recalibrationTables)
{
    int tempCOVarray[3];

    for(int i=0; i<reportTable.getNumRows(); i++)
    {
        tempCOVarray[0] = 0;
        tempCOVarray[1] = stoi(reportTable.get(i, RecalUtils::QUALITY_SCORE_COLUMN_NAME));

        string covName = reportTable.get(i, RecalUtils::COVARIATE_NAME_COLUMN_NAME);
        string covValue = reportTable.get(i, RecalUtils::COVARIATE_VALUE_COLUMN_NAME);
        Covariate * cov = covariates->getCovariateByParsedName(covName);
        assert(cov != nullptr);
        tempCOVarray[2] = cov->keyFromValue(covValue);

        (*recalibrationTables.getTableForCovariate(covName))[tempCOVarray[0]][tempCOVarray[1]][tempCOVarray[2]] = getRecalDatum(reportTable, i,
                                                                                                                             false);
    }
}

RecalDatum * RecalibrationReport::getRecalDatum(GATKReportTable &reportTable, int row, bool hasEstimatedQReportedColumn)
{
    long nObservations = stol(reportTable.get(row, RecalUtils::NUMBER_OBSERVATIONS_COLUMN_NAME));
    double nErrors = stod(reportTable.get(row, RecalUtils::NUMBER_ERRORS_COLUMN_NAME));

    // the estimatedQreported column only exists in the ReadGroup table
    double estimatedQReported = hasEstimatedQReportedColumn ? stod(reportTable.get(row, RecalUtils::ESTIMATED_Q_REPORTED_COLUMN_NAME)) :
            stoi(reportTable.get(row, RecalUtils::QUALITY_SCORE_COLUMN_NAME));

    RecalDatum * datum = new RecalDatum(nObservations, nErrors, 1);
    datum->setEstimatedQReported(estimatedQReported);
    return datum;
}

shared_ptr<RecalibrationTable> RecalibrationReport::getRecalibrationTables()
{
    return recalibrationTables;
}

shared_ptr<StandardCovariateList> RecalibrationReport::getCovariates()
{
    return covariates;
}

shared_ptr<QuantizationInfo> RecalibrationReport::getQuantizationInfo()
{
    return quantizationInfo;
}