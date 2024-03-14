//
// Created by 梦想家xixi on 2021/11/8.
//

#include <stdexcept>
#include "CigarOperator.h"

std::map<CigarOperator, CigarOperatorUtils> CigarOperatorUtils::cigarMap;

void CigarOperatorUtils::initial() {
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(M, CigarOperatorUtils(true, true, 'M')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(I, CigarOperatorUtils(true, false, 'I')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(D, CigarOperatorUtils(false, true, 'D')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(N, CigarOperatorUtils(false, true, 'N')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(S, CigarOperatorUtils(true, false, 'S')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(H, CigarOperatorUtils(false, false, 'H')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(P, CigarOperatorUtils(false, false, 'P')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(EQ, CigarOperatorUtils(true, true, '=')));
    cigarMap.insert(std::pair<CigarOperator, CigarOperatorUtils>(X, CigarOperatorUtils(true, true, 'X')));
}

bool CigarOperatorUtils::getConsumesReadBases(CigarOperator cigarOperator) {
    //return cigarMap.find(cigarOperator)->second.consumesReadBases;
    switch(cigarOperator) {
        case M: return true;
        case I: return true;
        case D: return false;
        case N: return false;
        case S: return true;
        case H: return false;
        case P: return false;
        case EQ: return true;
        case X: return true;
    }
}

bool CigarOperatorUtils::getConsumesReferenceBases(CigarOperator cigarOperator) {
    //return cigarMap.find(cigarOperator)->second.consumesReferenceBases;
    switch(cigarOperator) {
        case M: return true;
        case I: return false;
        case D: return true;
        case N: return true;
        case S: return false;
        case H: return false;
        case P: return false;
        case EQ: return true;
        case X: return true;
    }
}

CigarOperator CigarOperatorUtils::characterToEnum(int b) {
    switch(b) {
        case 61:
            return EQ;
        case 62:
        case 63:
        case 64:
        case 65:
        case 66:
        case 67:
        case 69:
        case 70:
        case 71:
        case 74:
        case 75:
        case 76:
        case 79:
        case 81:
        case 82:
        case 84:
        case 85:
        case 86:
        case 87:
        default:
            throw std::invalid_argument("Unrecognized CigarOperator" );
        case 68:
            return D;
        case 72:
            return H;
        case 73:
            return I;
        case 77:
            return M;
        case 78:
            return N;
        case 80:
            return P;
        case 83:
            return S;
        case 88:
            return X;
    }
}

CigarOperator CigarOperatorUtils::binaryToEnum(int i) {
    switch(i) {
        case 0:
            return M;
        case 1:
            return I;
        case 2:
            return D;
        case 3:
            return N;
        case 4:
            return S;
        case 5:
            return H;
        case 6:
            return P;
        case 7:
            return EQ;
        case 8:
            return X;
        default:
            throw std::invalid_argument("Unrecognized CigarOperator" );
    }
}

int CigarOperatorUtils::enumToBinary(CigarOperator e) {
    switch(e) {
        case M:
            return 0;
        case I:
            return 1;
        case D:
            return 2;
        case N:
            return 3;
        case S:
            return 4;
        case H:
            return 5;
        case P:
            return 6;
        case EQ:
            return 7;
        case X:
            return 8;
        default:
            throw std::invalid_argument("Unrecognized CigarOperator" );
    }
}

uint8_t CigarOperatorUtils::enumToCharacter(CigarOperator e) {
    return cigarMap.find(e)->second.character;
}

bool CigarOperatorUtils::isClipping(CigarOperator e) {
    return e == S || e == H;
}

bool CigarOperatorUtils::isIndel(CigarOperator e) {
    return e == I || e == D;
}

bool CigarOperatorUtils::isIndelOrSkippedRegion(CigarOperator e) {
    return e == N || isIndel(e);
}

bool CigarOperatorUtils::isAlignment(CigarOperator e) {
    return e == M || e == X || e == EQ;
}

bool CigarOperatorUtils::isPadding(CigarOperator e) {
    return e == P;
}
