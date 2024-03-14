//
// Created by 梦想家xixi on 2021/12/24.
//

#include "SAMTools_decode.h"
#include <iostream>
#include "SamFileHeaderMerger.h"

SAMFileHeader *SAMTools_decode::decode_samFileHeader(sam_hdr_t* header) {
    SAMFileHeader* samFileHeader = new SAMFileHeader();
    samFileHeader->init();
    setSequenceDictory(header, samFileHeader);
    setReadGroups(header, samFileHeader);
    setProgramRecords(header, samFileHeader);
    return samFileHeader;
}

void SAMTools_decode::setSequenceDictory(sam_hdr_t *header, SAMFileHeader * samFileHeader) {
    std::vector<SAMSequenceRecord> toAdd;
    int nref = sam_hdr_nref(header);
    for(int i = 0; i < nref; i++) {
        int length = (int)sam_hdr_tid2len(header, i);
        std::string name(sam_hdr_tid2name(header, i));
        SAMSequenceRecord samSequenceRecord(name, length);
        samSequenceRecord.setSequenceIndex(i);
        toAdd.emplace_back(samSequenceRecord);
    }
    samFileHeader->setSequenceDictionary(toAdd);
}

void SAMTools_decode::setReadGroups(sam_hdr_t * header, SAMFileHeader *samFileHeader) {
    std::vector<SAMReadGroupRecord> toAdd;
    int nRG = sam_hdr_count_lines(header, "RG");
    kstring_t tmp;
    ks_initialize(&tmp);
    for(int i = 0; i < nRG; i++) {
        sam_hdr_find_tag_pos(header, "RG", i, "ID", &tmp);
        std::string id(tmp.s);
        SAMReadGroupRecord samReadGroupRecord(id);
        ks_free(&tmp);
        setAttribute(header, "RG", i, "BC", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "CN", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "DS", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "DT", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "FO", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "KS", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "LB", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "PG", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "PI", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "PL", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "PM", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "PU", &tmp, &samReadGroupRecord);
        setAttribute(header, "RG", i, "SM", &tmp, &samReadGroupRecord);
        toAdd.emplace_back(samReadGroupRecord);
    }
    samFileHeader->setReadGroups(toAdd);
}

void SAMTools_decode::setAttribute(sam_hdr_t *h, const char *type, int pos, const char *key, kstring_t *ks,
                                   AbstractSAMHeaderRecord *abstractSamHeaderRecord) {
    ks_free(ks);
    ks_initialize(ks);
    if(sam_hdr_find_tag_pos(h, type, pos, key, ks) == 0) {
        std::string toAddKey(key);
        std::string value(ks->s);
        abstractSamHeaderRecord->setAttribute(toAddKey, value);
        ks_free(ks);
    }
    else {
        ks_free(ks);
        return;
    }
}

void SAMTools_decode::setProgramRecords(sam_hdr_t * header, SAMFileHeader *samFileHeader) {
    std::vector<SAMProgramRecord> toAdd;
    int nRG = sam_hdr_count_lines(header, "PG");
    kstring_t tmp;
    ks_initialize(&tmp);
    for(int i = 0; i < nRG; i++) {
        sam_hdr_find_tag_pos(header, "PG", i, "ID", &tmp);
        std::string id(tmp.s);
        SAMProgramRecord samProgramRecord(id);
        ks_free(&tmp);
        setAttribute(header, "PG", i, "PN", &tmp, &samProgramRecord);
        setAttribute(header, "PG", i, "CL", &tmp, &samProgramRecord);
        setAttribute(header, "PG", i, "PP", &tmp, &samProgramRecord);
        setAttribute(header, "PG", i, "DS", &tmp, &samProgramRecord);
        setAttribute(header, "PG", i, "VN", &tmp, &samProgramRecord);
        toAdd.emplace_back(samProgramRecord);
    }
    samFileHeader->setProgramRecords(toAdd);
}

SAMFileHeader *SAMTools_decode::merge_samFileHeaders(std::vector<sam_hdr_t *> headers) {
    std::vector<SAMFileHeader*> toDecode;
    for(sam_hdr_t* header : headers) {
        SAMFileHeader* samFileHeader = SAMTools_decode::decode_samFileHeader(header);
        toDecode.emplace_back(samFileHeader);
    }
    SamFileHeaderMerger merger = SamFileHeaderMerger();
    std::vector<SAMProgramRecord> pg = merger.mergeProgramGroups(toDecode);
    std::vector<SAMReadGroupRecord> rg = merger.mergeReadGroups(toDecode);
    SAMSequenceDictionary samSequenceDictionary = merger.mergeSequenceDictionaries(toDecode);
    SAMFileHeader* ret = new SAMFileHeader();
    ret->setProgramRecords(pg);
    ret->setSequenceDictionary(samSequenceDictionary.getSequences());
    ret->setReadGroups(rg);
    for(SAMFileHeader* samFileHeader1 : toDecode) {
        delete samFileHeader1;
    }
    return ret;
}



