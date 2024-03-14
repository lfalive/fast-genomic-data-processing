//
// Created by 梦想家xixi on 2021/11/8.
//

#ifndef MUTECT2CPP_MASTER_CIGAROPERATOR_H
#define MUTECT2CPP_MASTER_CIGAROPERATOR_H

#include <cstdint>
#include <map>

enum CigarOperator {M, I, D, N, S, H, P, X, EQ, _NULL}; // _NULL is added to represent null

class CigarOperatorUtils {
private:
    char character;
    bool consumesReadBases;
    bool consumesReferenceBases;

public:
    CigarOperatorUtils(bool consumesReadBases, bool consumesReferenceBases, char character) : consumesReadBases(consumesReadBases), consumesReferenceBases(consumesReferenceBases), character(character){}
    static std::map<CigarOperator, CigarOperatorUtils> cigarMap;
    static void initial();
    static const CigarOperator MATCH_OR_MISMATCH = M;
    static const CigarOperator INSERTION = I;
    static const CigarOperator DELETION = D;
    static const CigarOperator SKIPPED_REGION = N;
    static const CigarOperator SOFT_CLIP = S;
    static const CigarOperator HARD_CLIP = H;
    static const CigarOperator PADDING = P;
    static bool getConsumesReadBases(CigarOperator e);
    static bool getConsumesReferenceBases(CigarOperator e);
    static CigarOperator characterToEnum(int b);
    static CigarOperator binaryToEnum(int i);
    static int enumToBinary(CigarOperator e);
    static uint8_t enumToCharacter(CigarOperator e);
    static bool isClipping(CigarOperator e);
    static bool isIndel(CigarOperator e);
    static bool isIndelOrSkippedRegion(CigarOperator e);
    static bool isAlignment(CigarOperator e);
    static bool isPadding(CigarOperator e);
};


#endif //MUTECT2CPP_MASTER_CIGAROPERATOR_H
