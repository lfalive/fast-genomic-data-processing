//
// This class has all the static functionality for reading a recalibration report file into memory.
// Created by lhh on 4/1/22.
//

#ifndef MUTECT2CPP_MASTER_RECALIBRATIONREPORT_H
#define MUTECT2CPP_MASTER_RECALIBRATIONREPORT_H

#include "GATKReport.h"
#include "QuantizationInfo.h"
#include "RecalibrationArgumentCollection.h"
#include "RecalibrationTables.h"

class RecalibrationReport {
private:
    shared_ptr<GATKReportTable> argumentTable;
    shared_ptr<QuantizationInfo> quantizationInfo;
    shared_ptr<RecalibrationArgumentCollection> RAC;
    shared_ptr<StandardCovariateList> covariates;
    shared_ptr<RecalibrationTable> recalibrationTables; // quick access reference to tables

    static RecalibrationArgumentCollection* initializeArgumentCollectionTable(GATKReportTable & table);

    static RecalDatum* getRecalDatum(GATKReportTable& reportTable, int row, bool hasEstimatedQReportedColumn);

public:
    RecalibrationReport(char * recal_table);
    RecalibrationReport(GATKReport report);
    ~RecalibrationReport();

    shared_ptr<RecalibrationTable> getRecalibrationTables();

    shared_ptr<StandardCovariateList> getCovariates();

    shared_ptr<QuantizationInfo> getQuantizationInfo();

    /**
     * Parses the quantization table from the GATK Report and turns it into a map of original => quantized quality scores
     *
     * @param table the GATKReportTable containing the quantization mappings
     * @return an ArrayList with the quantization mappings from 0 to MAX_SAM_QUAL_SCORE
     */
    QuantizationInfo* initializeQuantizationTable(GATKReportTable & table);

    /**
     * Compiles the list of keys for the ReadGroup table and uses the shared parsing utility to produce the actual table
     *
     * @param reportTable            the GATKReport table containing data for this table
     * @param rgTable                the map representing this table
     */
    void parseReadGroupTable(GATKReportTable& reportTable, vector<RecalDatum *>& readGroupTable);

    /**
     *
     * Compiles the list of keys for the QualityScore table and uses the shared parsing utility to produce the actual table
     * @param reportTable            the GATKReport table containing data for this table
     * @param qualTable               the map representing this table
     */
    void parseQualityScoreTable(GATKReportTable& reportTable, TwoDimensionArray * qualTable);

    void parseAllCovariatesTable(GATKReportTable& reportTable, RecalibrationTable & recalibrationTables);

};


#endif //MUTECT2CPP_MASTER_RECALIBRATIONREPORT_H
