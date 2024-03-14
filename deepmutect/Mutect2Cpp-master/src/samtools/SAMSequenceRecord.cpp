//
// Created by 梦想家xixi on 2021/12/15.
//

#include "SAMSequenceRecord.h"


const std::string SAMSequenceRecord::SEQUENCE_NAME_TAG("SN");
const std::string SAMSequenceRecord::SEQUENCE_LENGTH_TAG("LN");
const std::string SAMSequenceRecord::MD5_TAG("M5");
const std::string SAMSequenceRecord::ASSEMBLY_TAG("AS");
const std::string SAMSequenceRecord::URI_TAG("UR");
const std::string SAMSequenceRecord::SPECIES_TAG("SP");
const std::string SAMSequenceRecord::DESCRIPTION_TAG("DS");
const std::string SAMSequenceRecord::RESERVED_RNEXT_SEQUENCE_NAME("=");
const std::regex SAMSequenceRecord::LEGAL_RNAME_PATTERN("[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*");


SAMSequenceRecord::SAMSequenceRecord(std::string &name, int sequenceLength) {
    mSequenceIndex = -1;
    mSequenceLength = 0;
    if(!name.empty()) {
        validateSequenceName(name);
        mSequenceName = name;
    }
    mSequenceLength = sequenceLength;
}

void SAMSequenceRecord::validateSequenceName(std::string &name) {
    if(!std::regex_match(name, LEGAL_RNAME_PATTERN))
        throw std::invalid_argument("Sequence name doesn't match regex");
}
