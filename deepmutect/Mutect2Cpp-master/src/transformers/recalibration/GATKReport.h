//
// Created by lhh on 4/1/22.
//

#ifndef MUTECT2CPP_MASTER_GATKREPORT_H
#define MUTECT2CPP_MASTER_GATKREPORT_H

#include <set>
#include <memory>
#include "GATKReportTable.h"

class GATKReport {
private:
    map<string, shared_ptr<GATKReportTable>> tables;

    /**
     * parse the Recalibration table file into a class
     */
    void parseReport(char * recal_table);

    /**
     * Adds a table, empty or populated, to the report
     *
     * @param table the table to add
     */
    void addTable(GATKReportTable* table);

public:
    GATKReport(char * recal_table);
    ~GATKReport();

    shared_ptr<GATKReportTable> getTable(std::string tableName);

    set<string> getReadGroups();
};


#endif //MUTECT2CPP_MASTER_GATKREPORT_H
