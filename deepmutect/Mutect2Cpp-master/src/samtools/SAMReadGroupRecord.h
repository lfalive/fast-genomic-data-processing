//
// Created by 梦想家xixi on 2021/12/22.
//

#ifndef MUTECT2CPP_MASTER_SAMREADGROUPRECORD_H
#define MUTECT2CPP_MASTER_SAMREADGROUPRECORD_H

#include <string>
#include <set>
#include "AbstractSAMHeaderRecord.h"
#include "HeaderRecordFactory.h"

class SAMReadGroupRecord : public AbstractSAMHeaderRecord, public HeaderRecordFactory{
private:
    std::string mReadGroupId;
public:
    static const std::string READ_GROUP_ID_TAG;
    static const std::string SEQUENCING_CENTER_TAG;
    static const std::string DESCRIPTION_TAG;
    static const std::string DATE_RUN_PRODUCED_TAG;
    static const std::string FLOW_ORDER_TAG;
    static const std::string KEY_SEQUENCE_TAG;
    static const std::string LIBRARY_TAG;
    static const std::string PROGRAM_GROUP_TAG;
    static const std::string PREDICTED_MEDIAN_INSERT_SIZE_TAG;
    static const std::string PLATFORM_TAG;
    static const std::string PLATFORM_MODEL_TAG;
    static const std::string PLATFORM_UNIT_TAG;
    static const std::string READ_GROUP_SAMPLE_TAG;
    static const std::string BARCODE_TAG;
    static const std::set<std::string> STANDARD_TAGS;
    explicit SAMReadGroupRecord(std::string & id);
    SAMReadGroupRecord(std::string & id, SAMReadGroupRecord & srcProgramRecord);
    ~SAMReadGroupRecord() override = default;
    std::string &  getReadGroupId();
    std::string & getId() override;
    AbstractSAMHeaderRecord* createRecord(std::string &newId, AbstractSAMHeaderRecord *record) override;
};


#endif //MUTECT2CPP_MASTER_SAMREADGROUPRECORD_H
