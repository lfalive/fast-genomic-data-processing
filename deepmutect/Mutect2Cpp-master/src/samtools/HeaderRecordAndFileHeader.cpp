//
// Created by 梦想家xixi on 2021/12/27.
//

#include "HeaderRecordAndFileHeader.h"

HeaderRecordAndFileHeader::HeaderRecordAndFileHeader(SAMFileHeader *samFileHeader,
                                                                          AbstractSAMHeaderRecord *headerRecord) : samFileHeader(samFileHeader), headerRecord(headerRecord) {}

SAMFileHeader * HeaderRecordAndFileHeader::getFileHeader() {
    return samFileHeader;
}

AbstractSAMHeaderRecord * HeaderRecordAndFileHeader::getHeaderRecord() {
    return headerRecord;
}
