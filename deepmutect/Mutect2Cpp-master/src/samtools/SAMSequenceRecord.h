//
// Created by 梦想家xixi on 2021/12/15.
//

#ifndef MUTECT2CPP_MASTER_SAMSEQUENCERECORD_H
#define MUTECT2CPP_MASTER_SAMSEQUENCERECORD_H

#include <string>
#include <regex>

class SAMSequenceRecord {
private:
    int mSequenceIndex;
    int mSequenceLength;
    std::string mSequenceName;
    static const std::string SEQUENCE_NAME_TAG;
    static const std::string SEQUENCE_LENGTH_TAG;
    static const std::string MD5_TAG;
    static const std::string ASSEMBLY_TAG;
    static const std::string URI_TAG;
    static const std::string SPECIES_TAG;
    static const std::string DESCRIPTION_TAG;
    static const std::string RESERVED_RNEXT_SEQUENCE_NAME;
    const int UNKNOWN_SEQUENCE_LENGTH = 0;
    static const std::regex LEGAL_RNAME_PATTERN;

public:
    SAMSequenceRecord(std::string & name, int sequenceLength);
    static void validateSequenceName(std::string & name);
    std::string & getSequenceName() {return mSequenceName;}
    void setSequenceIndex(int value) {mSequenceIndex = value;}
    int getSequenceLength() const {return mSequenceLength;}
    int getSequenceIndex() const {return mSequenceIndex;}
};


#endif //MUTECT2CPP_MASTER_SAMSEQUENCERECORD_H
