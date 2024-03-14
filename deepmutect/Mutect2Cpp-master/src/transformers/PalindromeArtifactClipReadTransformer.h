//
// Created by 梦想家xixi on 2022/2/22.
//

#ifndef MUTECT2CPP_MASTER_PALINDROMEARTIFACTCLIPREADTRANSFORMER_H
#define MUTECT2CPP_MASTER_PALINDROMEARTIFACTCLIPREADTRANSFORMER_H

#include "ReferenceCache.h"
#include "samtools/SAMRecord.h"

class PalindromeArtifactClipReadTransformer {
private:
    static const double MIN_FRACTION_OF_MATCHING_BASES;

    std::shared_ptr<ReferenceCache> referenceDataSource;

    SAMFileHeader* header;

    int minPalindromeSize;

public:
    PalindromeArtifactClipReadTransformer(std::shared_ptr<ReferenceCache>  referenceDataSource, SAMFileHeader* header,int minPalindromeSize);

    std::shared_ptr<SAMRecord> apply(const std::shared_ptr<SAMRecord> & read);

    // overrite apply() without using SAMRecord class for performance
    bam1_t * apply(bam1_t * read, sam_hdr_t * hdr);
};


#endif //MUTECT2CPP_MASTER_PALINDROMEARTIFACTCLIPREADTRANSFORMER_H
