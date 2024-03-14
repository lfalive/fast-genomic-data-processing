/**
 * The implementation of class Covariate , ContextCovariate
 * */

#include <iostream>
#include <vector>
#include <memory>
#include <cassert>
#include "ContextCovariate.h"
#include "htslib/sam.h"
#include "SamRead.h"

using namespace std;

void Covariate::recordValues(bam1_t *read, sam_hdr_t *readsHeader, Key *keys)
{
    cout << "This is Covariate" << endl;
}

int Covariate::maximumKeyValue()
{
    return 0;
}

ContextCovariate::ContextCovariate(BaseArgument & baseArgument):mismatchesContextSize(baseArgument.MismatchesContextSize), lowQualTail(baseArgument.LowQualityTail) {
    //this->LENGTH_BITS = 4;

    this->mismatchesKeyMask = createMask(mismatchesContextSize);
}
ContextCovariate::ContextCovariate(RecalibrationArgumentCollection &RAC) : mismatchesContextSize(RAC.MISMATCHES_CONTEXT_SIZE), lowQualTail(RAC.LOW_QUAL_TAIL)
{
    assert(mismatchesContextSize > 0);
    assert(mismatchesContextSize <= MAX_DNA_CONTEXT);

    mismatchesKeyMask = createMask(mismatchesContextSize);
}


int ContextCovariate::createMask(int contextSize) {
    int mask =0;
    // create 2*contextSize worth of bits
    for (int i = 0; i < contextSize; i++) {
        mask = (mask << 2) | 3;
    }
    // shift 4 bits to mask out the bits used to encode the length
    return mask << LENGTH_BITS;
}

void ContextCovariate::recordValues(bam1_t *read, sam_hdr_t *readsHeader, Key * keys) {
    int ReadLength = read->core.l_qseq;

    //---get the bases first, the bases has been transformed
    uint8_t * strandedClippedBases = getStrandedClippedBytes(read);

    vector<int>* mismatchKeys = contextWith(strandedClippedBases, ReadLength);

    if (mismatchKeys->size() != ReadLength)
    {
        for (int i=0; i<ReadLength; i++)
            keys[i].ContextKey = 0;
    }
    else{
        bool negativeStrand = bam_is_rev(read);

        for (int i=0; i<ReadLength; i++)
        {
            int readOffset = getStrandedOffset(negativeStrand, i, ReadLength);
            keys[readOffset].ContextKey = mismatchKeys->operator[](i);
        }
    }

    delete mismatchKeys;
    delete[] strandedClippedBases;
}

int ContextCovariate::getStrandedOffset(bool isNegativeStrand, int offset, int readLength) {
    return isNegativeStrand ? (readLength - offset -1) : offset;
}

vector<int> * ContextCovariate::contextWith(uint8_t * bases, int readLength){
    //int keys = new int(readLength);
    vector<int>* keys = new vector<int>;

    // if bases is NULL, there is a memory leak
    if(!bases)
	    return keys;

    keys->reserve(readLength);
    int i=0;

    // the first contextSize-1 bases will not have enough previous context
    for (; i< this->mismatchesContextSize-1 && i<readLength; i++)
        keys->push_back(-1);

    if (readLength < this->mismatchesContextSize)
        return keys;

    int newBaseOffset = ((this->mismatchesContextSize - 1) << 1) + LENGTH_BITS;

    // get (and add) the key for the context starting at the first base
    int currentKey = keyFromContext(bases, 0, this->mismatchesContextSize);
    keys->push_back(currentKey);

    // if the first key was -1 then there was an N in the context; figure out how many more consecutive contexts it affects
    int currentNPenalty = 0;
    if (currentKey == -1) {
        currentKey = 0;
        currentNPenalty = this->mismatchesContextSize - 1;
        int offset = newBaseOffset;
        while (bases[currentNPenalty] != 15) {  //---15 represent base N
            int baseIndex = simpleBaseToBaseIndex(bases[currentNPenalty]);
            currentKey |= (baseIndex << offset);
            offset -= 2;
            currentNPenalty--;
        }
    }

    for (int currentIndex = this->mismatchesContextSize; currentIndex < readLength; currentIndex++) {
        int baseIndex = simpleBaseToBaseIndex(bases[currentIndex]);
        if (baseIndex == -1) { // ignore non-ACGT bases
            currentNPenalty = this->mismatchesContextSize;
            currentKey = 0; // reset the key
        } else {
            // push this base's contribution onto the key: shift everything 2 bits, mask out the non-context bits, and add the new base and the length in
            currentKey = (currentKey >> 2) & this->mismatchesKeyMask;
            currentKey |= (baseIndex << newBaseOffset);
            currentKey |= this->mismatchesContextSize;
        }

        if (currentNPenalty == 0) {
            keys->push_back(currentKey);
        } else {
            currentNPenalty--;
            keys->push_back(-1);
        }
    }

    return keys;
}

int ContextCovariate::keyFromContext(uint8_t *dna, int start, int end) {
    int key = end - start;
    int bitOffset = LENGTH_BITS;
    for (int i=start; i<end; i++)
    {
        int baseIndex = simpleBaseToBaseIndex(dna[i]);
        if (baseIndex == -1)
            return -1;
        key |= (baseIndex << bitOffset);
        bitOffset += 2;
    }
    return key;
}

int ContextCovariate::keyFromContext(const char *dna, int start, int end)
{
    int key = end - start;
    int bitOffset = LENGTH_BITS;
    for (int i=start; i<end; i++)
    {
        int baseIndex = simpleBaseToBaseIndex(dna[i]);
        if (baseIndex == -1)
            return -1;
        key |= (baseIndex << bitOffset);
        bitOffset += 2;
    }
    return key;
}

int ContextCovariate::simpleBaseToBaseIndex(uint8_t base)       //--- change the function to a table?
{
    switch (base) {
        case 1:
            return 0;
        case 2:
            return 1;
        case 4:
            return 2;
        case 8:
            return 3;
        default:
            return -1;
    }
}

int ContextCovariate::simpleBaseToBaseIndex(char base)
{
    switch (base) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

//---please pay attention, each base is represented by 8 bits, not 4 bits in returned array
uint8_t * ContextCovariate::getStrandedClippedBytes(bam1_t *read) {
    read = clipLowQualEnds(read);

    if (bam_is_rev(read))
        return simpleReverseComplement(read);
    else{
	    if(!read->core.l_qseq)  // when clipLowQualEnds return an empty read
        {
            bam_destroy1(read);
            return nullptr;
        }

        uint8_t * bases = new uint8_t[read->core.l_qseq];	// bases is not nullptr when ReadLength == 0
        uint8_t * initialBases = bam_get_seq(read);
        for (int i=0; i<read->core.l_qseq; i++)
            bases[i] = bam_seqi(initialBases, i);
        return bases;
    }
}

uint8_t* ContextCovariate::simpleReverseComplement(bam1_t *read){
    int readLength = read->core.l_qseq;
    uint8_t * rcbases = new uint8_t[readLength];    //---Don't forget to delete it
    uint8_t basei;
    uint8_t * initialBases = bam_get_seq(read);

    for (int i=0; i<readLength; i++)
    {
        basei = bam_seqi(initialBases, readLength -i -1);
        rcbases[i] = simpleComplement(basei);
    }

    /*
    //store new bases into bam1_t struct
    for (int i=0; i<readLength; i++)
    {
        bam_set_seqi(bam_get_seq(read), i, rcbases[i]);
    }

    delete rcbases;
    return bam_get_seq(read);
     */

    return rcbases;
}

uint8_t ContextCovariate::simpleComplement(uint8_t base){
    switch (base) {
        case 1:
            return 8;
        case 2:
            return 4;
        case 4:
            return 2;
        case 8:
            return 1;
        default:
            return base;
    }
}

bam1_t* ContextCovariate::clipLowQualEnds(bam1_t *read){
    if (read->core.l_qseq == 0)
        return read;

    int readLength = read->core.l_qseq;
    int leftClipIndex = 0;
    int rightClipIndex = readLength - 1;


    uint8_t * qualities = getBaseQualities(read);

    // check how far we can clip both sides
    while (rightClipIndex >= 0 && qualities[rightClipIndex] <= lowQualTail){
        rightClipIndex--;
    }
    while (leftClipIndex < readLength && qualities[leftClipIndex] <= lowQualTail){
        leftClipIndex++;
    }

    // if the entire read should be clipped, then return an empty read.
    if (leftClipIndex > rightClipIndex) {
        return bam_init1();     //---to be confirmed
    }

    if(rightClipIndex < readLength - 1)
    {
        applyWriteNs(read, rightClipIndex + 1, readLength - 1);
    }

    if (leftClipIndex > 0)
    {
        applyWriteNs(read, 0, leftClipIndex - 1);
    }

    return read;
}

void ContextCovariate::applyWriteNs(bam1_t *read, int start, int end)
{
    for (int i=start; i<=end; i++)
        bam_set_seqi(bam_get_seq(read), i, seq_nt16_table['N']);
}
/*
int ContextCovariate::maximumKeyValue()
{
    int key = mismatchesContextSize;
    int bitOffset = LENGTH_BITS;
    for (int i=0; i<mismatchesContextSize; i++)
    {
        key |= (3 << bitOffset);
        bitOffset += 2;
    }
    return key;
}*/

int ContextCovariate::maximumKeyValue()
{
    return ContextCovariate::maximumKey;
}

int ContextCovariate::keyFromValue(string value)
{
    return keyFromContext(value.c_str(), 0, value.length());
}

string ContextCovariate::getName()
{
    return "Context";
}

string ContextCovariate::formatKey(int key)
{
    assert(key >= 0);
    int length = key & LENGTH_MASK; // the first bits represent the length (in bp) of the context
    int mask = 48; // use the mask to pull out bases
    int offset = LENGTH_BITS;

    unique_ptr<char[]> dna(new char[length+1]);
    for (int i = 0; i < length; i++) {
        int baseIndex = (key & mask) >> offset;
        dna[i] = baseIndexToSimpleBase(baseIndex);
        mask <<= 2; // move the mask over to the next 2 bits
        offset += 2;
    }
    dna[length] = '\0';
    return string(dna.get());
}

char ContextCovariate::baseIndexToSimpleBase(int baseIndex)
{
    switch (baseIndex) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            return '.';
    }
}

void ContextCovariate::printReads(bam1_t *read)
{
    int ReadLength = read->core.l_qseq;
    uint8_t * base = bam_get_seq(read);

    char RealBase;
    for(int i=0; i<ReadLength; i++)
    {

        char basei = bam_seqi(base, i);
        switch (basei) {
            case 1:
                RealBase = 'A';
                break;
            case 2:
                RealBase = 'C';
                break;
            case 4:
                RealBase = 'G';
                break;
            case 8:
                RealBase = 'T';
                break;
            case 15:
                RealBase = 'N';
                break;
            default:
                RealBase = '-';
                break;
        }
        cout << RealBase;
    }
    cout << endl;
}


