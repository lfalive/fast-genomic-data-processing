//
// Created by 梦想家xixi on 2021/12/22.
//

#ifndef MUTECT2CPP_MASTER_SAMFILEHEADER_H
#define MUTECT2CPP_MASTER_SAMFILEHEADER_H

#include "AbstractSAMHeaderRecord.h"
#include "SAMReadGroupRecord.h"
#include "SAMProgramRecord.h"
#include "SAMSequenceDictionary.h"
#include <set>
#include <vector>
#include <map>

enum SortOrder{
    SAMFileHeader_unsorted,
    SAMFileHeader_queryname,
    SAMFileHeader_coordinate,
    SAMFileHeader_duplicate,
    SAMFileHeader_unknown,
    SAMFileHeader_SortOrder_null
};

enum GroupOrder{
    SAMFileHeader_none,
    SAMFileHeader_query,
    SAMFileHeader_reference,
    SAMFileHeader_GroupOrder_null,
};

class SAMFileHeader : public AbstractSAMHeaderRecord{
public:
    static const std::string VERSION_TAG;
    static const std::string SORT_ORDER_TAG;
    static const std::string GROUP_ORDER_TAG;
    static const std::string CURRENT_VERSION;
    static const std::set<std::string> ACCEPTABLE_VERSIONS;
    static const std::set<std::string> STANDARD_TAGS;
    //需要init
    SAMFileHeader();
    void init();
    void setAttribute(std::string& key, std::string& value) override;
    void setTextHeader(char* text);
    int getSequenceIndex(std::string &basicString);
    void setSequenceDictionary(std::vector<SAMSequenceRecord> & toAdd);
    void setReadGroups(std::vector<SAMReadGroupRecord> & readGroups);
    void setProgramRecords(std::vector<SAMProgramRecord> & programRecords);
    std::vector<SAMProgramRecord> & getProgramRecords() {return mProgramRecords;}
    std::vector<SAMReadGroupRecord> & getReadGroupRecord() {return mReadGroups;}
    SAMSequenceDictionary& getSequenceDictionary() {return mSequenceDictionary;};

private:
    std::vector<SAMReadGroupRecord> mReadGroups;
    std::vector<SAMProgramRecord> mProgramRecords;
    std::map<std::string, SAMReadGroupRecord &> mReadGroupMap;
    std::map<std::string, SAMProgramRecord &> mProgramRecordMap;
    SAMSequenceDictionary mSequenceDictionary;
    std::vector<std::string> mComments;
    std::string textHeader;
    SortOrder sortOrder;
    GroupOrder groupOrder;
};


#endif //MUTECT2CPP_MASTER_SAMFILEHEADER_H
