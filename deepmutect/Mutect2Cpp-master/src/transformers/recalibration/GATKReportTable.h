/**
 * A helper class for output
 */

#ifndef GATK_REPORT_TABLE_H
#define GATK_REPORT_TABLE_H

#include <string>
#include <vector>
#include <map>
#include "parallel_hashmap/phmap.h"
#include "GATKReportColumn.h"
#include "SamRead.h"

using namespace std;

enum Sorting {
    SORT_BY_ROW,
    SORT_BY_COLUMN,
    DO_NOT_SORT
};



class GATKReportTable
{
private:

    string tableName;
    string tableDescription;
    Sorting sortingWay;

    vector<vector<string>> underlyingData;
    phmap::flat_hash_map<string, int> columnNameToIndex;
    vector<GATKReportColumn>* columnInfo;
    map<int, int> rowIdToIndex;

protected:
    enum TableDataHeaderFields {
        COLS = 2,
        ROWS = 3,
        FORMAT_START = 4
    };

    enum TableNameHeaderFields{
        NAME = 2,
        DESCRIPTION = 3
    };

public:
    const static string GATKREPORT_HEADER_PREFIX;
    const static string GATKREPORT_VERSION;
    const static int TablesSize = 5;    // in GATK, there should be a GATKReport class object containing 5 tables

    const static string GATKTABLE_HEADER_PREFIX;
    const static string SEPARATOR;
    const static string ENDLINE;
    const static string READGROUP_REPORT_TABLE_TITLE;
    const static string QUALITY_SCORE_REPORT_TABLE_TITLE;
    const static string ALL_COVARIATES_REPORT_TABLE_TITLE;

    GATKReportTable(string tableName, string tableDescription, int numColumns, Sorting sortingWay);

    GATKReportTable(ifstream & reader);

    ~GATKReportTable();

    string getTableName();

    int getNumRows();

    string & get(int rowIndex, int columnIndex);

    string & get(int rowIndex, string columnName);

    /**
     * Add a column to the report and the format string used to display the data
     * @param columnName
     * @param format
     *//*
    void addColumn(string columnName, string format);

    *//**
     * Set the value for a given position in the table.
     *//*
    void set(int rowID, string columnName, int value);

    void set(int rowIndex, int colIndex, int value);

    void expandTo(int rowIndex, bool updateRowIdMap);

    *//**
     * Get the number of columns in this table
     *//*
     int getNumColumns();*/

    void static outputArgumentTable(ofstream & output, BaseArgument & baseArgument);

    void static inputArgumentTable(ifstream & input, BaseArgument & baseArgument);
};

#endif