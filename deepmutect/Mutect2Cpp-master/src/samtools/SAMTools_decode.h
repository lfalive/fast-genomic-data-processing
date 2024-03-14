//
// Created by 梦想家xixi on 2021/12/24.
//

#ifndef MUTECT2CPP_MASTER_SAMTOOLS_DECODE_H
#define MUTECT2CPP_MASTER_SAMTOOLS_DECODE_H

#include "SAMFileHeader.h"
#include "htslib/sam.h"

class SAMTools_decode {
public:
    static SAMFileHeader* decode_samFileHeader(sam_hdr_t* header);

    // Don't forget to free the returned SAMFileHeader*
    static SAMFileHeader* merge_samFileHeaders(std::vector<sam_hdr_t *> headers);

private:
    static void setSequenceDictory(sam_hdr_t * header, SAMFileHeader* samFileHeader);
    static void setReadGroups(sam_hdr_t * header, SAMFileHeader* samFileHeader);
    static void setAttribute(sam_hdr_t *h, const char *type, int pos, const char *key, kstring_t *ks, AbstractSAMHeaderRecord* abstractSamHeaderRecord);
    static void setProgramRecords(sam_hdr_t *  header, SAMFileHeader* samFileHeader);
};


#endif //MUTECT2CPP_MASTER_SAMTOOLS_DECODE_H
