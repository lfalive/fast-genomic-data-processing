//
// Created by 梦想家xixi on 2021/11/15.
//

#ifndef MUTECT2CPP_MASTER_DANGLINGCHAINMERGEHELPER_H
#define MUTECT2CPP_MASTER_DANGLINGCHAINMERGEHELPER_H

#include "MultiDeBruijnVertex.h"
#include "Cigar.h"
#include <utility>
#include <vector>
#include <utility>

class DanglingChainMergeHelper {
public:
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> danglingPath;
    std::vector<std::shared_ptr<MultiDeBruijnVertex>> referencePath;
    std::shared_ptr<uint8_t[]> danglingPathString;
    std::shared_ptr<uint8_t[]> referencePathString;
    int danglingPathStringLength;
    int referencePathStringLength;
    std::shared_ptr<Cigar> cigar;

	DanglingChainMergeHelper(std::vector<std::shared_ptr<MultiDeBruijnVertex>> danglingPath, std::vector<std::shared_ptr<MultiDeBruijnVertex>> referencePath, std::shared_ptr<uint8_t[]> danglingPathString, int danglingPathStringLength,
                             std::shared_ptr<uint8_t[]> referencePathString,  int referencePathStringLength, std::shared_ptr<Cigar>  cigar) : danglingPath(std::move(danglingPath)), referencePath(std::move(referencePath)), danglingPathString(std::move(danglingPathString)),
                             referencePathString(std::move(referencePathString)), danglingPathStringLength(danglingPathStringLength), referencePathStringLength(referencePathStringLength), cigar(std::move(cigar)) {}
};


#endif //MUTECT2CPP_MASTER_DANGLINGCHAINMERGEHELPER_H
