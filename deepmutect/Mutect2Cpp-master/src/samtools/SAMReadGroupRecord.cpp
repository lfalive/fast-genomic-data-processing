//
// Created by 梦想家xixi on 2021/12/22.
//

#include "SAMReadGroupRecord.h"

const std::string SAMReadGroupRecord::READ_GROUP_ID_TAG = "ID";
const std::string SAMReadGroupRecord::SEQUENCING_CENTER_TAG = "CN";
const std::string SAMReadGroupRecord::DESCRIPTION_TAG = "DS";
const std::string SAMReadGroupRecord::DATE_RUN_PRODUCED_TAG = "DT";
const std::string SAMReadGroupRecord::FLOW_ORDER_TAG = "FO";
const std::string SAMReadGroupRecord::KEY_SEQUENCE_TAG = "KS";
const std::string SAMReadGroupRecord::LIBRARY_TAG = "LB";
const std::string SAMReadGroupRecord::PROGRAM_GROUP_TAG = "PG";
const std::string SAMReadGroupRecord::PREDICTED_MEDIAN_INSERT_SIZE_TAG = "PI";
const std::string SAMReadGroupRecord::PLATFORM_TAG = "PL";
const std::string SAMReadGroupRecord::PLATFORM_MODEL_TAG = "PM";
const std::string SAMReadGroupRecord::PLATFORM_UNIT_TAG = "PU";
const std::string SAMReadGroupRecord::READ_GROUP_SAMPLE_TAG = "SM";
const std::string SAMReadGroupRecord::BARCODE_TAG = "BC";
const std::set<std::string> SAMReadGroupRecord::STANDARD_TAGS{"ID", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", "PL", "PM", "PU", "SM", "BC"};

SAMReadGroupRecord::SAMReadGroupRecord(std::string &id) : mReadGroupId(id){
}

SAMReadGroupRecord::SAMReadGroupRecord(std::string &id, SAMReadGroupRecord &srcProgramRecord) : mReadGroupId(id){
    for(std::pair<std::string, std::string> tmp : srcProgramRecord.getAttributes()) {
        this->setAttribute(tmp.first, tmp.second);
    }
}

std::string &SAMReadGroupRecord::getReadGroupId() {
    return mReadGroupId;
}

std::string &SAMReadGroupRecord::getId() {
    return mReadGroupId;
}

AbstractSAMHeaderRecord *SAMReadGroupRecord::createRecord(std::string &newId, AbstractSAMHeaderRecord *record) {
    return new SAMReadGroupRecord(newId, *(SAMReadGroupRecord*)record);
}
