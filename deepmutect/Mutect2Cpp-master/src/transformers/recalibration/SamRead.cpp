/**
 * The implementation of SamRecord.h
 */

#include <iostream>
#include "SamRead.h"
#include "string.h"
//#include "RecalUtils.h"

#define MAPPING_QUALITY_UNAVAILABLE 255
#define NO_ALIGNMENT_REFERENCE_NAME "*"
#define NO_ALIGNMENT_START 0
#define UNSET_POSITION 0

const char* GetReferenceName(const sam_hdr_t *h, bam1_t *read)
{
    return sam_hdr_tid2name(h, read->core.tid);
}

hts_pos_t GetStart(bam1_t *read)
{
    return read->core.pos;
}

//---how to get the end of a read?
int GetEnd(bam1_t *read)
{
    return bam_endpos(read) - 1;
}

bool IsUnmapped(const sam_hdr_t *h, bam1_t *read)
{
    const char * refName = GetReferenceName(h, read);
    return (read->core.flag & BAM_FUNMAP) != 0 ||
           strcmp(refName, NO_ALIGNMENT_REFERENCE_NAME) == 0 ||
           GetStart(read) == NO_ALIGNMENT_START;
}

int getFragmentLength(bam1_t *read)
{
    return read->core.isize;
}

bool isPaired(bam1_t *read)
{
    return (read->core.flag & BAM_FPAIRED) != 0;
}

bool hasWellDefinedFragmentSize(const sam_hdr_t *h, bam1_t *read)
{
    if (getFragmentLength(read) == 0)
        return false;
    if (!isPaired(read))
        return false;
    if (IsUnmapped(h, read) || mateIsUnmapped(h, read))
        return false;
    if (isReverseStrand(read) == mateIsReverseStrand(read))
        return false;
    if (isReverseStrand(read))
        return GetEnd(read) > getMateStart(h, read);
    else
        return GetStart(read) <= getMateStart(h, read) + getFragmentLength(read);
}

bool isReverseStrand(bam1_t *read)
{
    return (read->core.flag & BAM_FREVERSE) != 0;
}

bool mateIsReverseStrand(bam1_t *read)
{
    return (read->core.flag & BAM_FMREVERSE) != 0;
}

bool getMateUnmappedFlag(bam1_t *read)
{
    if (!getReadPairedFlag(read))
        throw "Inappropriate call if not paired read";
    return getMateUnmappedFlagUnchecked(read);
}

bool getMateUnmappedFlagUnchecked(bam1_t *read)
{
    return (read->core.flag & BAM_FMUNMAP) != 0;
}

bool mateIsUnmapped(const sam_hdr_t *h, bam1_t *read)
{
    if (!isPaired(read))
        throw "Cannot get mate information for an unpaired read";
    return getMateUnmappedFlag(read) || getMateReferenceName(h, read) == NULL ||
            strcmp(getMateReferenceName(h, read), NO_ALIGNMENT_REFERENCE_NAME) == 0 ||
            getMateAlignmentStart(read) == NO_ALIGNMENT_START;
}

const char* getMateReferenceName(const sam_hdr_t *h, bam1_t *read)
{
    return sam_hdr_tid2name(h, read->core.mtid);
}

bool getReadPairedFlag(bam1_t * read)
{
    return (read->core.flag & BAM_FPAIRED) != 0;
}

bool getProperPairFlagUnchecked(bam1_t * read)
{
    return (read->core.flag & BAM_FPROPER_PAIR) != 0;
}

int getMateStart(const sam_hdr_t *h, bam1_t *read)
{
    if(mateIsUnmapped(h, read))
        return UNSET_POSITION;
    return getMateAlignmentStart(read);
}

int getMateAlignmentStart(bam1_t * read)
{
    return read->core.mpos;
}

bool isInsideRead(bam1_t *read, int referenceCoordinate)
{
    return referenceCoordinate >= GetStart(read) && referenceCoordinate <= GetEnd(read);
}

uint8_t * getBaseQualities(bam1_t * read)
{
    return bam_get_qual(read);
}

uint8_t getBaseQuality(bam1_t * read, int i)
{
    uint8_t * qualities = bam_get_qual(read);
    return qualities[i];
}

bool getSecondOfPairFlagUnchecked(bam1_t * read)
{
    return (read->core.flag & BAM_FREAD2) != 0;
}

int getSoftStart(bam1_t * read)
{
    int softStart = GetStart(read);
    uint32_t* cigar = bam_get_cigar(read);
    for (unsigned i=0; i<read->core.n_cigar; i++)
    {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CSOFT_CLIP)
            softStart -= bam_cigar_oplen(cigar[i]);
        else if (op != BAM_CHARD_CLIP)
            break;
    }
    return softStart;
}

int getSoftEnd(bam1_t * read)
{
    bool foundAlignedBase = false;
    int softEnd = GetEnd(read);
    uint32_t* cigar = bam_get_cigar(read);
    for (unsigned i=read->core.n_cigar-1; i>=0; i--)
    {
        int op = bam_cigar_op(cigar[i]);

        if (op == BAM_CSOFT_CLIP)
            softEnd += bam_cigar_oplen(cigar[i]);
        else if (op != BAM_CHARD_CLIP){
            foundAlignedBase = true;
            break;
        }
    }
    if (!foundAlignedBase)
        softEnd = GetEnd(read);

    return softEnd;
}

bool consumesReadBases(uint32_t cigarElement)
{
    return bam_cigar_type(bam_cigar_op(cigarElement)) & 1;
}

bool consumesReferenceBases(uint32_t cigarElement)
{
    return bam_cigar_type(bam_cigar_op(cigarElement)) & 2;
}

void PrintBases(bam1_t * read)
{
    uint8_t * bases = bam_get_seq(read);
    char base;
    for (int i = 0; i < read->core.l_qseq; ++i) {
        base = bam_seqi(bases, i);
        switch (base) {
            case 1:
                std::cout << 'A';
                break;
            case 2:
                std::cout << 'C';
                break;
            case 4:
                std::cout << 'G';
                break;
            case 8:
                std::cout << 'T';
                break;
            case 15:
                std::cout << 'N';
                break;
            default:
                std::cout << '-';
                break;
        }

    }
    std::cout << std::endl;
}

void PrintQualities(bam1_t * read)
{
    uint8_t * qualities = bam_get_qual(read);
    for (int i = 0; i < read->core.l_qseq; ++i) {
        std::cout << (int)qualities[i] << " ";
    }
    std::cout << std::endl;
}

void PrintCigar(bam1_t * read)
{
    uint32_t * cigar = bam_get_cigar(read);
    for (unsigned i = 0; i < read->core.n_cigar; ++i) {
        std::cout << bam_cigar_oplen(cigar[i]) << " ";

        switch (bam_cigar_op(cigar[i])) {
            case BAM_CMATCH:
                std::cout << "M" << "\t";
                break;
            case BAM_CINS:
                std::cout << "I" << "\t";
                break;
            case BAM_CDEL:
                std::cout << "D" << "\t";
                break;
            case BAM_CREF_SKIP:
                std::cout << "N" << "\t";
                break;
            case BAM_CSOFT_CLIP:
                std::cout << "S" << "\t";
                break;
            case BAM_CHARD_CLIP:
                std::cout << "H" << "\t";
                break;
        }

    }
    std::cout << std::endl;


}

void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                             hts_pos_t *rlen, hts_pos_t *qlen)
{
    int k;
    *rlen = *qlen = 0;
    for (k = 0; k < n_cigar; ++k) {
        int type = bam_cigar_type(bam_cigar_op(cigar[k]));
        int len = bam_cigar_oplen(cigar[k]);
        if (type & 1) *qlen += len;
        if (type & 2) *rlen += len;	
    }
}
