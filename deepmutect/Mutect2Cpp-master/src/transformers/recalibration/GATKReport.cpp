//
// Created by lhh on 4/1/22.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include "GATKReport.h"
#include "RecalUtils.h"

GATKReport::GATKReport(char *recal_table)
{
    parseReport(recal_table);
}

GATKReport::~GATKReport()
{
}

void GATKReport::parseReport(char *recal_table)
{
    assert(recal_table != nullptr);
    std::ifstream input(recal_table);
    std::string GATKReportHeader;
    if(input.is_open())
    {
        try {
            getline(input, GATKReportHeader);
            if(GATKReportHeader.empty())
                throw (std::string(recal_table) + "is empty !");
        } catch (const std::string msg){
            std::cerr<< msg << std::endl;
        }

        size_t found = GATKReportHeader.find_last_of(":");

        int nTables = stoi(GATKReportHeader.substr(found + 1));
        for(int i=0; i<nTables; i++)
        {
            addTable(new GATKReportTable(input));
        }

    } else {
        throw "cannot open the recalibration table file";
    }
}

void GATKReport::addTable(GATKReportTable * table)
{
    tables.insert(std::pair<string, shared_ptr<GATKReportTable>>(table->getTableName(), shared_ptr<GATKReportTable>(table)));
}

shared_ptr<GATKReportTable> GATKReport::getTable(std::string tableName)
{
    return tables[tableName];
}

set<string> GATKReport::getReadGroups()
{
    GATKReportTable* reportTable = getTable(RecalUtils::READGROUP_REPORT_TABLE_TITLE).get();
    set<string> readGroups;
    for(int i=0; i<reportTable->getNumRows(); i++)
    {
        readGroups.insert(reportTable->get(i, 0));
    }
    return readGroups;
}
