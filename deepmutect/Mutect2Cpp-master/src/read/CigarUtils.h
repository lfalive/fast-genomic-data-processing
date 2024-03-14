//
// Created by 梦想家xixi on 2021/11/25.
//

#ifndef MUTECT2CPP_MASTER_CIGARUTILS_H
#define MUTECT2CPP_MASTER_CIGARUTILS_H

#include "cigar/Cigar.h"
#include "SWParameters.h"
#include "SmithWatermanAlignment.h"

class CigarUtils {
public:
    static std::shared_ptr<Cigar> calculateCigar(const std::shared_ptr<uint8_t[]>& refSeq, int refLength, const std::shared_ptr<uint8_t[]>& altSeq, int altLength);
    static const SWParameters NEW_SW_PARAMETERS;

    // In Mutect2 and HaplotypeCaller reads are realigned to their *best* haplotypes, which is very different from a generic alignment.
    // The {@code NEW_SW_PARAMETERS} penalize a substitution error more than an indel up to a length of 9 bases!
    // Suppose, for example, that a read has a single substitution error, say C -> T, on its last base.  Those parameters
    // would prefer to extend a deletion until the next T on the reference is found in order to avoid the substitution, which is absurd.
    // Since these parameters are for aligning a read to the biological sequence we believe it comes from, the parameters
    // we choose should correspond to sequencer error.  They *do not* have anything to do with the prevalence of true variation!
    static const SWParameters ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS;
    static std::shared_ptr<Cigar> leftAlignCigarSequentially(std::shared_ptr<Cigar> & cigar, const std::shared_ptr<uint8_t[]>& refSeq, int refLength, const std::shared_ptr<uint8_t[]>& readSeq, int readLength, int refIndex, int readIndex);
    static bool isGood(const std::shared_ptr<Cigar> & c);
    static bool isGood(bam1_t * read);
    static bool containsNOperator(std::vector<CigarElement> cigarElements);
    static bool containsNOperator(uint32_t n_cigar, uint32_t* cigarArray);

private:
    static const int SW_PAD = 10;
    static bool isSWFailure(SmithWatermanAlignment* alignment);
    static bool hasConsecutiveIndels(std::vector<CigarElement> & elems);
    static bool hasConsecutiveIndels(uint32_t * cigarArray, uint32_t n_cigar);
    static bool startsWithDeletionIgnoringClips(const std::vector<CigarElement> & elems);
    static bool startsWithDeletionIgnoringClips(uint32_t * cigarArray, uint32_t n_cigar);
};


#endif //MUTECT2CPP_MASTER_CIGARUTILS_H
