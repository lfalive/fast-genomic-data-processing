//
// Created by 梦想家xixi on 2021/12/16.
//

#ifndef MUTECT2CPP_MASTER_SAMSEQUENCEDICTIONARY_H
#define MUTECT2CPP_MASTER_SAMSEQUENCEDICTIONARY_H

#include "SAMSequenceRecord.h"
#include "parallel_hashmap/phmap.h"
#include <map>
#include <vector>
#include <list>

class SAMSequenceDictionary {
private:
    std::vector<SAMSequenceRecord> mSequences;
    phmap::flat_hash_map<std::string, SAMSequenceRecord> mSequenceMap;

public:
    SAMSequenceDictionary() = default;
    explicit SAMSequenceDictionary(std::list<SAMSequenceRecord> & toAdd);
    void addSequence(SAMSequenceRecord sequenceRecord);
    SAMSequenceRecord & getSequence(const std::string& name);
    int getSequenceIndex(std::string & sequenceName);
    std::vector<SAMSequenceRecord> & getSequences() {return mSequences;}
};


#endif //MUTECT2CPP_MASTER_SAMSEQUENCEDICTIONARY_H
