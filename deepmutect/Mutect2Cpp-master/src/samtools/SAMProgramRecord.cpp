//
// Created by 梦想家xixi on 2021/12/22.
//

#include "SAMProgramRecord.h"

const std::string SAMProgramRecord::PROGRAM_GROUP_ID_TAG = "ID";
const std::string SAMProgramRecord::PROGRAM_NAME_TAG = "PN";
const std::string SAMProgramRecord::PROGRAM_VERSION_TAG = "VN";
const std::string SAMProgramRecord::COMMAND_LINE_TAG = "CL";
const std::string SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG = "PP";

SAMProgramRecord::SAMProgramRecord(std::string & programGroupId) : mProgramGroupId(programGroupId){
}

SAMProgramRecord::SAMProgramRecord(std::string &id, SAMProgramRecord &srcProgramRecord) : mProgramGroupId(id){
    for(std::pair<std::string, std::string> tmp : srcProgramRecord.getAttributes()) {
        this->setAttribute(tmp.first, tmp.second);
    }
}

std::string &SAMProgramRecord::getProgramGroupId() {
    return mProgramGroupId;
}

std::string &SAMProgramRecord::getId() {
    return mProgramGroupId;
}

AbstractSAMHeaderRecord *SAMProgramRecord::createRecord(std::string &newId, AbstractSAMHeaderRecord *record) {
    return new SAMProgramRecord(newId, *(SAMProgramRecord*)record);
}
