/**
 * column information within a GATK report table
 */

#ifndef GATK_REPORT_COLUMN_H
#define GATK_REPORT_COLUMN_H

#include <string>
using namespace std;
class GATKReportColumn{
private:
    string columnName;
    string format;
    string dataType;    //---in Java, type of this variable is GATKReportDataType
    int maxWidth = 0;

public:
    GATKReportColumn(string columnName, string format);

    string fromFormatString(string format);
};


#endif