/**
 * The implementation of GATKReportTable class
 */

#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "GATKReportTable.h"


const string GATKReportTable::GATKREPORT_HEADER_PREFIX = "#:GATKReport.";
const string GATKReportTable::GATKREPORT_VERSION = "v1.1";
const string GATKReportTable::SEPARATOR = ":";
const string GATKReportTable::GATKTABLE_HEADER_PREFIX = "#:GATKTable";
const string GATKReportTable::ENDLINE = ":;";
const string GATKReportTable::READGROUP_REPORT_TABLE_TITLE = "RecalTable0";
const string GATKReportTable::QUALITY_SCORE_REPORT_TABLE_TITLE = "RecalTable1";
const string GATKReportTable::ALL_COVARIATES_REPORT_TABLE_TITLE = "RecalTable2";

GATKReportTable::GATKReportTable(string tableName, string tableDescription, int numColumns, Sorting sortingWay)
{
    this->tableName = tableName;
    this->tableDescription = tableDescription;
    this->sortingWay = sortingWay;

    //columnInfo = new vector<GATKReportColumn>;


}

inline void split_str(string const & str, const char delim, vector<string> & out)
{
    stringstream s(str);

    string s2;
    while(getline(s, s2, delim))
    {
        if(!s2.empty())
            out.push_back(s2);
    }

}

GATKReportTable::GATKReportTable(ifstream &reader)
{
    vector<string> tableData, tableNameData;
    string line;

    getline(reader, line);
    split_str(line, ':', tableData);

    getline(reader, line);
    split_str(line, ':', tableNameData);

    // parse the header fields
    this->tableName = tableNameData[NAME];
    assert(tableNameData.size() >= DESCRIPTION);
    this->tableDescription = tableNameData.size() <= DESCRIPTION ? "" : tableNameData[DESCRIPTION];

    // initialize the data
    int nColumns = stoi(tableData[COLS]);
    int nRows = stoi(tableData[ROWS]);
    underlyingData.reserve(nRows);
    columnInfo = new vector<GATKReportColumn>;
    columnInfo->reserve(nColumns);
    columnNameToIndex.reserve(nColumns);

    // when reading from a file, the row ID mapping is just the index
    for(int i=0; i<nRows; i++)
    {
        rowIdToIndex.insert(std::pair<int, int>(i, i));
    }

    // read the column names
    string columnLine;
    getline(reader, columnLine);
    vector<string> columnNames;
    split_str(columnLine, '\t', columnNames);

    // Put in columns using the format string from the header
    for(int i=0; i<nColumns; i++)
    {
        columnNameToIndex.insert(std::pair<string, int>(columnNames[i], i));
    }

    // fill in the table
    for(int i=0; i<nRows; i++)
    {
        string dataLine;
        getline(reader, dataLine);
        vector<string> lineSplits;
        split_str(dataLine, '\t', lineSplits);
        assert(lineSplits.size() == nColumns);

        underlyingData.emplace_back(lineSplits);
    }

    // read the line between two tables
    getline(reader, line);
}



GATKReportTable::~GATKReportTable()
{
    delete columnInfo;
}

string GATKReportTable::getTableName()
{
    return tableName;
}

int GATKReportTable::getNumRows()
{
    return underlyingData.size();
}

string &GATKReportTable::get(int rowIndex, int columnIndex)
{
    assert(rowIndex >= 0 && rowIndex < underlyingData.size());
    assert(columnIndex >= 0 && columnIndex < underlyingData[0].size());
    return underlyingData[rowIndex][columnIndex];
}

string & GATKReportTable::get(int rowIndex, string columnName)
{
    return get(rowIndex, columnNameToIndex[columnName]);
}

/*void GATKReportTable::addColumn(string columnName, string format)
{
    columnNameToIndex.insert(std::make_pair(columnName, columnInfo->size()));
    columnInfo->push_back(GATKReportColumn(columnName, format));
}

void GATKReportTable::set(int rowID, string columnName, int value)
{
    if (rowIdToIndex.count(rowID) == 0){
        rowIdToIndex.insert(std::make_pair(rowID, underlyingData.size()));
        expandTo(underlyingData.size(), false);
    }
    set(rowIdToIndex.operator[](rowID), columnNameToIndex.operator[](columnName), value);
}

void GATKReportTable::set(int rowIndex, int colIndex, int value)
{
    expandTo(rowIndex, true);
    assert(rowIndex >=0 && colIndex >= 0 && colIndex >= getNumColumns());
    GATKReportColumn column = columnInfo->operator[](colIndex);
    //---TODO: finish the GATKReportTable class 2021.4.8

}*/

/*void GATKReportTable::expandTo(int rowIndex, bool updateRowIdMap)
{
    int currentSize = underlyingData.size();
    if (rowIndex >= currentSize)
    {
        int numNewRows = rowIndex - currentSize + 1;
        for (int i=0; i<numNewRows; i++)
        {
            if (updateRowIdMap)
                rowIdToIndex.insert(make_pair(currentSize, currentSize));
            underlyingData.push_back(vector<int>(getNumColumns()));
        }
    }
}*/

/*
int GATKReportTable::getNumColumns()
{
    return columnInfo->size();
}
*/

//TODO: make output and input more elegant
void GATKReportTable::outputArgumentTable(ofstream & output, BaseArgument & baseArgument) //TODO: replace some value to the parameter entered
{
    output << GATKReportTable::GATKTABLE_HEADER_PREFIX << ":" << 2 << ":" << 17 << ":" << "%s:%s" << GATKReportTable::ENDLINE << endl;
    output << GATKReportTable::GATKTABLE_HEADER_PREFIX << GATKReportTable::SEPARATOR << "Arguments" << GATKReportTable::SEPARATOR << "Recalibration argument collection values used in this run" << endl;
    output << "Argument" << "\t\t\t" << "Value" << endl;
    output << "binary_tag_name" << "\t\t\t" << "null" << endl;  // an optional parameter
    output << "covariate" << "\t\t\t" << "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate" << endl;
    output << "default-platform" << "\t\t" << "null" << endl;  // a hidden advanced parameter
    output << "deletions_default_quality" << "\t" << 45 << endl;  // an optional parameter
    output << "default_platform" << "\t\t" << "null" << endl;   // a hidden advanced parameter
    output << "indels_context_size" << "\t\t" << 3 << endl; // an optional parameter
    output << "insertions_default_quality" << "\t" << 45 << endl;
    output << "low_quality_tail" << "\t\t" << baseArgument.LowQualityTail << endl;
    output << "maximum_cycle_value" << "\t\t" << baseArgument.MaximumCycleValue << endl;
    output << "mismatches_context_size" << "\t\t" << baseArgument.MismatchesContextSize << endl;
    output << "mismatches_default_quality" << "\t" << baseArgument.MismatchesDefaultQuality << endl;
    output << "no_standard_covs" << "\t\t" << false << endl;
    output << "quantizing_levels" << "\t\t" << baseArgument.QuantizingLevels << endl;
    output << "recalibration_report" << "\t\t" << "null" << endl;
    output << "run_without_dbsnp" << "\t\t" << false << endl;
    output << "solid_nocall_strategy" << "\t\t" << "THROW_EXCEPTION" << endl;
    output << "solid_recal_mode" << "\t\t" << "SET_Q_ZERO" << endl;
    output << endl;
}

void GATKReportTable::inputArgumentTable(ifstream &input, BaseArgument &baseArgument)
{
    string ArgumentTableHeader;
    string TableDescription;
    input >> ArgumentTableHeader;
    input >> TableDescription;



}