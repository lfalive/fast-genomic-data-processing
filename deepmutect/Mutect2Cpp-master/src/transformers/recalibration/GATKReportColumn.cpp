/**
 * The implementation of the GATKReportColumn class
 */

#include "GATKReportColumn.h"

GATKReportColumn::GATKReportColumn(string columnName, string format)
{
    this->columnName = columnName;
    this->maxWidth = columnName.size();
    if (format == "")
    {
        this->format = "%s";
        this->dataType = "Unknown";
    }
    else
    {
        this->format = format;
        this->dataType = fromFormatString(format);  //
    }
}

string GATKReportColumn::fromFormatString(string format)
{
    if (format == "%d" || format == "%D")
        return "Integer";
    else if (format == "%b")
        return "Boolean";
    else if (format == "%c")
        return "Character";
    else if (format == "%f")
        return "Decimal";
    else if (format == "%s")
        return "String";
    else
        return "NULL";
}