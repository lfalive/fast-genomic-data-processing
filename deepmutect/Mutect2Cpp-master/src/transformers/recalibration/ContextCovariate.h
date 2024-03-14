/**
 *  Created by lhh, 2021.3.3
 * */

#ifndef CONTEXT_COVARIATE_H
#define CONTEXT_COVARIATE_H

#include <vector>
#include <string>
#include <map>
#include "htslib/sam.h"
#include "SamRead.h"
#include "RecalibrationArgumentCollection.h"
using namespace std;

struct Key{
    int ContextKey;
    int CycleKey;
};

class Covariate{
public:
    virtual ~Covariate() = default;   // virtual destructor is needed

    virtual void recordValues(bam1_t *read, sam_hdr_t *readsHeader, Key * keys);    //---what's the type of keys?
    virtual int maximumKeyValue();
    virtual int keyFromValue(string value) = 0;
    virtual string getName() = 0;
    //virtual string formatKey(int key) = 0;
};

class ContextCovariate: public Covariate{
private:
    int mismatchesContextSize;
    uint8_t lowQualTail;   //---in Java, its type is byte
    static const int LENGTH_BITS = 4;
    static const int LENGTH_MASK = 15;
    static const int MAX_DNA_CONTEXT = 13;
    int mismatchesKeyMask;


    /**
     * print a read alignment to stdout, for debugging
     * @param read
     */
    void printReads(bam1_t *read);
public:
    const static int maximumKey = 1011;    //in GATK, this value is calculated using maximumKeyValue()

    ContextCovariate(BaseArgument & baseArgument);

    ContextCovariate(RecalibrationArgumentCollection & RAC);

    static int createMask(int contextSize);

    void recordValues(bam1_t *read, sam_hdr_t *readsHeader, Key * keys);
    uint8_t* getStrandedClippedBytes(bam1_t *read);

    bam1_t* clipLowQualEnds(bam1_t *read);

    /**
     * transform the bases of [start, end] to N
     */
    void applyWriteNs(bam1_t *read, int start, int end);


    uint8_t* simpleReverseComplement(bam1_t *read);

    /**
     * Return the complement (A <-> T or C <-> G) of a base, or the specified base if it can't be complemented (i.e. an ambiguous base).
     *
     * @param base the base [AaCcGgTt]
     * @return the complementary base, or the input base if it's not one of the understood ones
     */
    uint8_t simpleComplement(uint8_t base);

    //int * contextWith(uint8_t * bases, int readLength);
    vector<int> * contextWith(uint8_t * bases, int readLength);

    /**
    * Creates a int representation of a given dna string.
    *
    * @param dna    the dna sequence. Pay attention! the base must be 4 bit encoded(1 2 4 or 8)
    * @param start  the start position in the byte array (inclusive)
    * @param end    the end position in the array (exclusive)
    * @return the key representing the dna sequence
    */
    static int keyFromContext(uint8_t * dna, int start, int end);

    static int keyFromContext(const char * dna, int start, int end);

    // each base is encoded as 1, 2, 4 or 8
    static int simpleBaseToBaseIndex(uint8_t base);

    // each base is represented as a char ('A' 'C' 'T' 'G')
    static int simpleBaseToBaseIndex(char base);

    /**
    * Helper method: computes the correct offset to use in computations of covariate values.
    * @param isNegativeStrand is the read on the negative strand
    * @param offset 0-based index of the base in the read
    * @param readLength length of the read
    * @return
    */
    static int getStrandedOffset(bool isNegativeStrand, int offset, int readLength);

    int maximumKeyValue();

    int keyFromValue(string value);

    string getName();

    static string formatKey(int key);

    /**
     * Converts a base index to a simple base
     *
     * @param baseIndex 0, 1, 2, 3
     * @return A, C, G, T, or '.' if the index can't be understood
     */
    static char baseIndexToSimpleBase(int baseIndex);
};



#endif
