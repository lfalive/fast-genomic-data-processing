//
// Created by 梦想家xixi on 2021/12/27.
//

#ifndef MUTECT2CPP_MASTER_HEADERRECORDANDFILEHEADER_H
#define MUTECT2CPP_MASTER_HEADERRECORDANDFILEHEADER_H

#include "SAMFileHeader.h"

class HeaderRecordAndFileHeader {
private:
    SAMFileHeader * samFileHeader;
    AbstractSAMHeaderRecord * headerRecord;

public:
    HeaderRecordAndFileHeader(SAMFileHeader * samFileHeader, AbstractSAMHeaderRecord * headerRecord);
    SAMFileHeader * getFileHeader();
    AbstractSAMHeaderRecord * getHeaderRecord();

};


#endif //MUTECT2CPP_MASTER_HEADERRECORDANDFILEHEADER_H
