/**
 * A helper class to cache the reference sequence
 */

#ifndef REFERENCE_CACHE_H
#define REFERENCE_CACHE_H

#include "htslib/faidx.h"
#include "samtools/SAMFileHeader.h"

class ReferenceCache
{
private:
    char * bases;
    int tid;
    hts_pos_t start;
    hts_pos_t end;
    hts_pos_t len{};
    SAMFileHeader* header;
    faidx_t * fai;

public:
    ReferenceCache(char * refName, SAMFileHeader* header, int tid);

    ~ReferenceCache();

    void setTid(int tid);

    // get the name of cached chromosome
    std::string getContig();

    // delete all the elements in the cache
    void clear();

    void advanceLoad();

    /**
     * Get a single base from the reference cache
     */
    char getBase(hts_pos_t pos);

    std::shared_ptr<uint8_t[]> getSubsequenceAt(int tid, int start, int stop, int & length);
};


#endif
