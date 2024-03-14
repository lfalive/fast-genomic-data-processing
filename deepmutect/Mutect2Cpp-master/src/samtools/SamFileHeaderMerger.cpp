//
// Created by 梦想家xixi on 2021/12/27.
//

#include "SamFileHeaderMerger.h"
#include <vector>
#include <list>

char SamFileHeaderMerger::INT_TO_BASE36[36]{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};



bool SamFileHeaderMerger::mergeHeaderRecords(std::vector<HeaderRecordAndFileHeader> &headerRecords, HeaderRecordFactory* headerRecordFactory,
                                             std::set<std::string> &idsThatAreAlreadyTaken,
                                             std::map<SAMFileHeader*, std::map<std::string, std::string>> &idTranslationTable,
                                             std::vector<AbstractSAMHeaderRecord*> &result) {
    std::map<std::string, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute>> idToRecord;
    for(HeaderRecordAndFileHeader& pair : headerRecords) {
        AbstractSAMHeaderRecord* record = pair.getHeaderRecord();
        SAMFileHeader* header = pair.getFileHeader();
        std::string recordId = record->getId();
        if(idToRecord.find(recordId) == idToRecord.end()) {
            idToRecord.insert(std::pair<std::string, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute>>(recordId, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute>()));
        }
        std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute>& recordsWithSameId = idToRecord.at(recordId);

        if(recordsWithSameId.find(record) == recordsWithSameId.end()) {
            recordsWithSameId.insert(std::pair<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>>(record, std::vector<SAMFileHeader*>()));
        }
        std::vector<SAMFileHeader*>& fileHeaders = recordsWithSameId.at(record);
        fileHeaders.emplace_back(header);
    }
    bool hasCollisions = false;
    std::map<std::string, std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute>>::iterator miter;
    for(miter = idToRecord.begin(); miter != idToRecord.end(); miter++) {
        std::string recordId = miter->first;
        std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute> & recordsWithSameId = miter->second;
        std::map<AbstractSAMHeaderRecord*, std::vector<SAMFileHeader*>, compare_Attribute>::iterator iter;
        for(iter = recordsWithSameId.begin(); iter != recordsWithSameId.end(); iter++) {
            AbstractSAMHeaderRecord* record = iter->first;
            std::vector<SAMFileHeader*> & fileHeaders = iter->second;
            std::string newId;
            if(idsThatAreAlreadyTaken.find(recordId) == idsThatAreAlreadyTaken.end()) {
                newId = recordId;
                idsThatAreAlreadyTaken.insert(recordId);
                ++recordCounter;
            } else {
                hasCollisions = true;
                while(idsThatAreAlreadyTaken.find(newId = recordId + "." + positiveFourDigitBase36Str(recordCounter++)) != idsThatAreAlreadyTaken.end());
                idsThatAreAlreadyTaken.insert(newId);
            }
            for(SAMFileHeader* fileHeader : fileHeaders) {
                if(idTranslationTable.find(fileHeader) == idTranslationTable.end()) {
                    idTranslationTable.insert(std::pair<SAMFileHeader*, std::map<std::string, std::string>>(fileHeader, std::map<std::string, std::string>()));
                }
                std::map<std::string, std::string> & readerTranslationTable = idTranslationTable.at(fileHeader);
                readerTranslationTable.insert(std::pair<std::string, std::string>(recordId, newId));
            }
            AbstractSAMHeaderRecord* toAdd = headerRecordFactory->createRecord(newId, record);
            toDelete.emplace_back(toAdd);
            result.emplace_back(toAdd);
        }
    }
    return hasCollisions;
}

std::string SamFileHeaderMerger::positiveFourDigitBase36Str(int leftOver) {
    if(leftOver == 0) {
        return "0";
    } else {
        std::string ret;
        while (leftOver > 0) {
            int valueIndex = leftOver % 36;
            ret += INT_TO_BASE36[valueIndex];
            leftOver /= 36;
        }
        std::reverse(ret.begin(), ret.end());
        return ret;
    }
}

std::vector<SAMProgramRecord> SamFileHeaderMerger::mergeProgramGroups(std::vector<SAMFileHeader *> headers) {
    std::vector<SAMProgramRecord> overallResult;
    std::set<std::string> idsThatAreAlreadyTaken;
    std::vector<HeaderRecordAndFileHeader> programGroupsLeftToProcess;
    for(SAMFileHeader* header : headers) {
        for(SAMProgramRecord& programRecord : header->getProgramRecords()) {
            if(idsThatAreAlreadyTaken.find(programRecord.getId()) != idsThatAreAlreadyTaken.end()) {
                throw std::invalid_argument("Input file contains more than one PG with the same id");
            } else {
                idsThatAreAlreadyTaken.insert(programRecord.getId());
            }
            programGroupsLeftToProcess.emplace_back(HeaderRecordAndFileHeader(header, &programRecord));
        }
        idsThatAreAlreadyTaken.clear();
    }
    recordCounter = 0;
    std::vector<HeaderRecordAndFileHeader> currentProgramGroups;
    for(std::vector<HeaderRecordAndFileHeader>::iterator programGroupsLeftToProcessIterator = programGroupsLeftToProcess.begin();
    programGroupsLeftToProcessIterator != programGroupsLeftToProcess.end(); ) {
        HeaderRecordAndFileHeader pair = *programGroupsLeftToProcessIterator;
        if(pair.getHeaderRecord()->getAttribute(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG).empty()) {
            programGroupsLeftToProcess.erase(programGroupsLeftToProcessIterator);
            currentProgramGroups.emplace_back(pair);
        } else {
            programGroupsLeftToProcessIterator++;
        }
    }
    while(!currentProgramGroups.empty()) {
        std::vector<AbstractSAMHeaderRecord*> currentResult;
        std::string tmp = "tmp";
        HeaderRecordFactory* factory = new SAMProgramRecord(tmp);
        hasProgramGroupCollisions |= mergeHeaderRecords(currentProgramGroups, factory, idsThatAreAlreadyTaken, samProgramGroupIdTranslation, currentResult);
        delete (SAMProgramRecord*)factory;
        for(AbstractSAMHeaderRecord* samHeaderRecord : currentResult) {
            overallResult.emplace_back(*(SAMProgramRecord*)samHeaderRecord);
        }
        currentProgramGroups = translateIds(currentProgramGroups, samProgramGroupIdTranslation, toDelete, false);
        programGroupsLeftToProcess = translateIds(programGroupsLeftToProcess, samProgramGroupIdTranslation, toDelete, true);

        std::vector<HeaderRecordAndFileHeader> programGroupsToProcessNext;
        for(std::vector<HeaderRecordAndFileHeader>::iterator programGroupsLeftToProcessIterator = programGroupsLeftToProcess.begin();
        programGroupsLeftToProcessIterator != programGroupsLeftToProcess.end();) {
            bool flag = true;
            std::string ppIdOfRecordLeftToProcess = programGroupsLeftToProcessIterator->getHeaderRecord()->getAttribute(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG);
            for(HeaderRecordAndFileHeader justProcessedPair : currentProgramGroups) {
                std::string idJustProcessed = justProcessedPair.getHeaderRecord()->getId();
                if(programGroupsLeftToProcessIterator->getFileHeader() == justProcessedPair.getFileHeader() && ppIdOfRecordLeftToProcess == idJustProcessed) {
                    programGroupsToProcessNext.emplace_back(*programGroupsLeftToProcessIterator);
                    programGroupsLeftToProcessIterator = programGroupsLeftToProcess.erase(programGroupsLeftToProcessIterator);
                    flag = false;
                    break;
                }
            }
            if(flag) {
                programGroupsLeftToProcessIterator++;
            }
        }
        currentProgramGroups = programGroupsToProcessNext;
    }
    if(!programGroupsLeftToProcess.empty()) {
        throw std::invalid_argument("program groups weren't processed. Do their PP ids point to existing PGs?");
    }
    std::sort(overallResult.begin(), overallResult.end(), compareById);
    return overallResult;
}

std::vector<HeaderRecordAndFileHeader> SamFileHeaderMerger::translateIds(const std::vector<HeaderRecordAndFileHeader>& programGroups,
                                  std::map<SAMFileHeader *, std::map<std::string, std::string>> idTranslationTable, std::vector<AbstractSAMHeaderRecord*> &todelete, bool translatePpIds) {
    std::vector<HeaderRecordAndFileHeader> result;
    for(HeaderRecordAndFileHeader pair : programGroups) {
        SAMProgramRecord& record = *(SAMProgramRecord*)pair.getHeaderRecord();
        std::string id = record.getId();
        std::string ppId = record.getAttribute(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG);
        SAMFileHeader* header = pair.getFileHeader();
        std::map<std::string, std::string> translations = idTranslationTable.find(header) != idTranslationTable.end() ? idTranslationTable.at(header) : std::map<std::string, std::string>();
        SAMProgramRecord* translatedRecord = nullptr;
        if(!translations.empty()) {
            std::string translatedId = translations.find(id) != translations.end() ? translations.at(id) : "";
            std::string translatedPpId = translatePpIds ? (translations.find(ppId) != translations.end() ? translations.at(ppId) : "") : "";
            bool needToTranslateId = !translatedId.empty() && translatedId != id;
            bool needToTranslatePpId = !translatedPpId.empty() && translatedPpId != ppId;
            if(needToTranslateId && needToTranslatePpId) {
                translatedRecord = new SAMProgramRecord(translatedId, record);
                translatedRecord->setAttribute(
                        const_cast<std::string &>(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG), translatedPpId);
            } else if (needToTranslateId) {
                translatedRecord = new SAMProgramRecord(translatedId, record);
            } else if (needToTranslatePpId) {
                translatedRecord = new SAMProgramRecord(id, record);
                translatedRecord->setAttribute(
                        const_cast<std::string &>(SAMProgramRecord::PREVIOUS_PROGRAM_GROUP_ID_TAG), translatedPpId);
            }
        }
        if(translatedRecord != nullptr) {
            result.emplace_back(HeaderRecordAndFileHeader(header, translatedRecord));
            todelete.emplace_back(translatedRecord);
        } else {
            result.emplace_back(pair);
        }
    }
    return result;
}

SamFileHeaderMerger::SamFileHeaderMerger() {
    hasProgramGroupCollisions = false;
    hasReadGroupCollisions = false;
}

bool SamFileHeaderMerger::compareByAttributes(AbstractSAMHeaderRecord *a, AbstractSAMHeaderRecord *b){
    std::map <std::string, std::string> a_map = a->getAttributes();
    std::map <std::string, std::string> b_map = b->getAttributes();
    for(std::pair<std::string, std::string> aPair : a_map) {
        if(b_map.find(aPair.first) == b_map.end())
            return false;
        if(aPair.second != b_map.at(aPair.first)) {
            return aPair.second < b_map.at(aPair.first);
        }
    }
    return false;
}

SAMSequenceDictionary
SamFileHeaderMerger::mergeSequences(SAMSequenceDictionary mergeIntoDict, SAMSequenceDictionary mergeFromDict) {
    std::list<SAMSequenceRecord> holder;

    std::list<SAMSequenceRecord> resultingDict;

    for(const SAMSequenceRecord& sequenceRecord : mergeIntoDict.getSequences()) {
        resultingDict.emplace_back(sequenceRecord);
    }
    int prevloc = -1;
    SAMSequenceRecord* previouslyMerged = nullptr;

    for(SAMSequenceRecord& sequenceRecord : mergeFromDict.getSequences()) {
        int loc = getIndexOfSequenceName(resultingDict, sequenceRecord.getSequenceName());
        if(loc == -1) {
            holder.emplace_back(sequenceRecord);
        } else if (prevloc > loc) {
            throw std::invalid_argument("Cannot merge sequence dictionaries");
        } else {
            std::list<SAMSequenceRecord>::iterator iter = std::next(resultingDict.begin(), loc);
            resultingDict.insert((resultingDict.begin()++) , holder.begin(), holder.end());
            prevloc = loc + holder.size();
            previouslyMerged = &sequenceRecord;
            holder.clear();
        }
    }
    if(!holder.empty()) {
        for(SAMSequenceRecord& samProgramRecord : holder) {
            resultingDict.emplace_back(samProgramRecord);
        }
    }
    return SAMSequenceDictionary(resultingDict);
}

int SamFileHeaderMerger::getIndexOfSequenceName(std::vector<SAMSequenceRecord> &list, const std::string& sequenceName) {
    for(int i = 0; i < list.size(); ++i) {
        if(list[i].getSequenceName() == sequenceName) {
            return i;
        }
    }
    return -1;
}

int SamFileHeaderMerger::getIndexOfSequenceName(std::list<SAMSequenceRecord> &list, const std::string& sequenceName) {
    int count = 0;
    for(std::list<SAMSequenceRecord>::iterator iter = list.begin(); iter != list.end(); iter++) {
        if(iter->getSequenceName() == sequenceName) {
            return count;
        } else {
            count++;
        }
    }
    return -1;
}

void SamFileHeaderMerger::createSequenceMapping(std::vector<SAMFileHeader *> headers,
                                                SAMSequenceDictionary &masterDictionary) {
    std::list<std::string> resultingDictStr;
    for(SAMSequenceRecord r : masterDictionary.getSequences()) {
        resultingDictStr.emplace_back(r.getSequenceName());
    }
    for(SAMFileHeader* header : headers) {
        std::map<int, int> seqMap;
        SAMSequenceDictionary& dict = header->getSequenceDictionary();
        for(SAMSequenceRecord rec : dict.getSequences()) {
            seqMap.insert(std::pair<int, int>(rec.getSequenceIndex(), indexOf(resultingDictStr, rec.getSequenceName())));
        }
        samSeqDictionaryIdTranslationViaHeader.insert(std::pair<SAMFileHeader*, std::map<int, int>> (header, seqMap));
    }
}

int SamFileHeaderMerger::indexOf(std::list<std::string> &stringList, std::string &stringToFind) {
    int count = 0;
    std::list<std::string>::iterator iter = stringList.begin();
    while(iter != stringList.end()) {
        if(*iter == stringToFind) {
            break;
        }
        count++;
        iter++;
    }
    return count;
}

SAMSequenceDictionary SamFileHeaderMerger::mergeSequenceDictionaries(std::vector<SAMFileHeader *>& headers) {
    SAMSequenceDictionary sequences = SAMSequenceDictionary();
    for(SAMFileHeader* header : headers) {
        SAMSequenceDictionary& currentSequences = header->getSequenceDictionary();
        sequences = mergeSequences(sequences, currentSequences);
    }
    createSequenceMapping(headers, sequences);
    return sequences;
}

std::vector<SAMReadGroupRecord> SamFileHeaderMerger::mergeReadGroups(std::vector<SAMFileHeader *> headers) {
    std::set<std::string> idsThatAreAlreadyTaken;
    std::vector<HeaderRecordAndFileHeader> readGroupsToProcess;
    for(SAMFileHeader* header : headers) {
        for(SAMReadGroupRecord& readGroup : header->getReadGroupRecord()) {
            if(idsThatAreAlreadyTaken.find(readGroup.getId()) != idsThatAreAlreadyTaken.end()) {
                throw std::invalid_argument("Input file: contains more than one RG with the same id");
            } else {
                idsThatAreAlreadyTaken.insert(readGroup.getId());
            }
            readGroupsToProcess.emplace_back(HeaderRecordAndFileHeader(header, &readGroup));
        }
        idsThatAreAlreadyTaken.clear();
    }
    std::vector<AbstractSAMHeaderRecord*> result;

    recordCounter = 0;
    std::string tmp = "tmp";
    HeaderRecordFactory * factory = new SAMReadGroupRecord(tmp);
    hasReadGroupCollisions = mergeHeaderRecords(readGroupsToProcess, factory, idsThatAreAlreadyTaken, samProgramGroupIdTranslation, result);
    delete factory;
    std::vector<SAMReadGroupRecord> ret;
    for(AbstractSAMHeaderRecord* samHeaderRecord : result) {
        ret.emplace_back(*(SAMReadGroupRecord*)samHeaderRecord);
    }
    std::sort(ret.begin(), ret.end(), compareById);
    return ret;
}

SamFileHeaderMerger::~SamFileHeaderMerger() {
    for(AbstractSAMHeaderRecord* samHeaderRecord : toDelete) {
        delete samHeaderRecord;
    }
}
