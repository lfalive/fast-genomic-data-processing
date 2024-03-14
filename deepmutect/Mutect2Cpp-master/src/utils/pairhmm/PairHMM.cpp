//
// Created by lhh on 4/22/22.
//

#include <cassert>
#include "PairHMM.h"

void PairHMM::initialize(const std::vector<std::shared_ptr<Haplotype>> &haplotypes,
                         const std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>> &perSampleReadList, int readMaxLength,
                         int haplotypeMaxLength) {
    initialize(readMaxLength, haplotypeMaxLength);
}

void PairHMM::initialize(int readMaxLength, int haplotypeMaxLength)
{
    assert(haplotypeMaxLength > 0);

    maxHaplotypeLength = haplotypeMaxLength;
    maxReadLength = readMaxLength;

    // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
    paddedMaxReadLength = readMaxLength + 1;
    paddedMaxHaplotypeLength = haplotypeMaxLength + 1;

    constantsAreInitialized = false;
    initialized = true;
}

PairHMM::~PairHMM() {

}
