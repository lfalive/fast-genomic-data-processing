//
// Created by 梦想家xixi on 2021/12/16.
//

#include "SAMSequenceDictionary.h"

void SAMSequenceDictionary::addSequence(SAMSequenceRecord sequenceRecord) {
    if(mSequenceMap.find(sequenceRecord.getSequenceName()) != mSequenceMap.end()) {
        throw std::invalid_argument("Cannot add sequence that already exists in SAMSequenceDictionary");
    } else {
        sequenceRecord.setSequenceIndex(mSequences.size());
        mSequences.emplace_back(sequenceRecord);
        mSequenceMap.insert(std::pair<std::string, SAMSequenceRecord>(sequenceRecord.getSequenceName(), sequenceRecord));
    }
}

SAMSequenceRecord & SAMSequenceDictionary::getSequence(const std::string& name) {
    if(mSequenceMap.find(name) == mSequenceMap.end())
        throw std::invalid_argument("not found in provided dictionary");
    return mSequenceMap.at(name);
}

int SAMSequenceDictionary::getSequenceIndex(std::string &sequenceName) {
    try {
        SAMSequenceRecord record = mSequenceMap.at(sequenceName);
        return record.getSequenceIndex();
    } catch(const std::out_of_range &e) {
        return -1;
    }
}

SAMSequenceDictionary::SAMSequenceDictionary(std::list<SAMSequenceRecord> &toAdd) {
    for(SAMSequenceRecord samSequenceRecord : toAdd) {
        addSequence(samSequenceRecord);
    }
}
