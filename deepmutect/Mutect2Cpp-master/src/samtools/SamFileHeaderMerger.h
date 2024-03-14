//
// Created by 梦想家xixi on 2021/12/27.
//

#ifndef MUTECT2CPP_MASTER_SAMFILEHEADERMERGER_H
#define MUTECT2CPP_MASTER_SAMFILEHEADERMERGER_H

#include "SAMFileHeader.h"
#include "HeaderRecordAndFileHeader.h"
#include <list>
struct compare_Attribute{
    bool operator() (const AbstractSAMHeaderRecord * a, const AbstractSAMHeaderRecord * b) const {
        std::map <std::string, std::string> a_map = a->getAttributesNochange();
        std::map <std::string, std::string> b_map = b->getAttributesNochange();
        for(std::pair<std::string, std::string> aPair : a_map) {
            if(b_map.find(aPair.first) == b_map.end())
                return false;
            if(aPair.second != b_map.at(aPair.first)) {
                return aPair.second < b_map.at(aPair.first);
            }
        }
        return false;
    }
};

class SamFileHeaderMerger {
private:
    int recordCounter;
    static char INT_TO_BASE36[36];
    bool hasProgramGroupCollisions;
    bool hasReadGroupCollisions;
    std::vector<AbstractSAMHeaderRecord*> toDelete;
    std::map<SAMFileHeader*, std::map<std::string, std::string>> samProgramGroupIdTranslation;
    std::map<SAMFileHeader*, std::map<int, int>> samSeqDictionaryIdTranslationViaHeader;

    bool mergeHeaderRecords(std::vector<HeaderRecordAndFileHeader> &headerRecords, HeaderRecordFactory* headerRecordFactory, std::set<std::string> & idsThatAreAlreadyTaken, std::map<SAMFileHeader*, std::map<std::string, std::string>> & idTranslationTable,
                            std::vector<AbstractSAMHeaderRecord*> & result);



    std::vector<HeaderRecordAndFileHeader> translateIds(const std::vector<HeaderRecordAndFileHeader>& programGroups, std::map<SAMFileHeader*, std::map<std::string, std::string>>, std::vector<AbstractSAMHeaderRecord*> &todelete, bool translatePpIds);

    static bool compareById (AbstractSAMHeaderRecord& a, AbstractSAMHeaderRecord& b) {return a.getId() < b.getId();}

    static bool compareByAttributes(AbstractSAMHeaderRecord* a, AbstractSAMHeaderRecord* b) ;

    SAMSequenceDictionary mergeSequences(SAMSequenceDictionary mergeIntoDict, SAMSequenceDictionary mergeFromDict);

    static int getIndexOfSequenceName(std::vector<SAMSequenceRecord> &list, const std::string& sequenceName);

    static int getIndexOfSequenceName(std::list<SAMSequenceRecord> &list, const std::string& sequenceName);

    void createSequenceMapping(std::vector<SAMFileHeader*> headers, SAMSequenceDictionary & masterDictionary);

    int indexOf(std::list<std::string> & stringList, std::string & stringToFind);

public:
    std::string positiveFourDigitBase36Str(int leftOver);

    SamFileHeaderMerger();

    ~SamFileHeaderMerger();

    std::vector<SAMProgramRecord> mergeProgramGroups(std::vector<SAMFileHeader*> headers);

    SAMSequenceDictionary mergeSequenceDictionaries(std::vector<SAMFileHeader*>& headers);

    std::vector<SAMReadGroupRecord> mergeReadGroups(std::vector<SAMFileHeader*> headers);
};


#endif //MUTECT2CPP_MASTER_SAMFILEHEADERMERGER_H
