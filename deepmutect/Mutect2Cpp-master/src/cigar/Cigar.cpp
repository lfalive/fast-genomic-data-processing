//
// Created by 梦想家xixi on 2021/11/9.
//

#include <stdexcept>
#include "Cigar.h"
#include "string"
#include "ReadUtils.h"

Cigar::Cigar(std::vector<CigarElement>& cigarElements) {
    for(CigarElement e : cigarElements)
        this->cigarElements.emplace_back(e);
}

std::vector<CigarElement> & Cigar::getCigarElements() {
    return cigarElements;
}

CigarElement &Cigar::getCigarElement(int i) {
    return cigarElements[i];
}

void Cigar::add(CigarElement cigarElement) {
    cigarElements.emplace_back(cigarElement);
}

int Cigar::getReferenceLength() {
    int length = 0;
    for(CigarElement e : cigarElements) {
        switch (e.getOperator()) {
            case M:
            case D:
            case N:
            case EQ:
            case X:
                length += e.getLength();
        }
    }
    return length;
}

int Cigar::getPaddedReferenceLength() {
    int length = 0;
    for(CigarElement e : cigarElements) {
        switch (e.getOperator()) {
            case M:
            case D:
            case N:
            case EQ:
            case X:
            case P:
                length += e.getLength();
        }
    }
    return length;
}

int Cigar::getReadLength(const std::vector<CigarElement> & cigarElements) {
    int length = 0;
    for(CigarElement e : cigarElements) {
        if(CigarOperatorUtils::getConsumesReadBases(e.getOperator())) {
            length += e.getLength();
        }
    }

    return length;
}

int Cigar::getReadLength(uint32_t n_cigar, uint32_t *cigarArray)
{
    int length = 0;
    for(uint32_t i=0; i<n_cigar; i++)
    {
        if(ReadUtils::consumesReadBases(cigarArray[i])){
            length += bam_cigar_oplen(cigarArray[i]);
        }
    }
    return length;
}

int Cigar::getReadLength() {
    return getReadLength(cigarElements);
}

Cigar* Cigar::fromCigarOperators(std::vector<CigarElement> &cigarElements) {
    if(cigarElements.empty())
        throw std::invalid_argument("cigarOperators is null");
    else {
        std::vector<CigarElement> cigarElementList;
        int j;
        for(int i = 0; i < cigarElements.size(); i = j) {
            CigarOperator currentOp = cigarElements[i].getOperator();
            for(j = i + 1; j < cigarElements.size() && cigarElements[j].getOperator() == currentOp; ++j){}
            cigarElementList.emplace_back(CigarElement(j - i, currentOp));
        }
        Cigar* res = new Cigar(cigarElementList);
        return res;
    }
}

bool Cigar::isRealOperator(CigarOperator op) {
    return op == M || op == EQ || op == X || op == I || op == D || op == N;
}

bool Cigar::isRealOperator(uint32_t op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CINS || op == BAM_CDEL || op == BAM_CREF_SKIP;
}

bool Cigar::isInDelOperator(CigarOperator op) {
    return CigarOperatorUtils::isIndel(op);
}
bool Cigar::isInDelOperator(uint32_t op) {
    return op == BAM_CINS || op == BAM_CDEL;
}

bool Cigar::isClippingOperator(CigarOperator op) {
    return CigarOperatorUtils::isClipping(op);
}

bool Cigar::isPaddingOperator(CigarOperator op) {
    return CigarOperatorUtils::isPadding(op);
}

bool Cigar::isPaddingOperator(uint32_t op)
{
    return op == BAM_CPAD;
}

bool Cigar::isLeftClipped() {
    return !isEmpty() && isClippingOperator(cigarElements[0].getOperator());
}

bool Cigar::isRightClipped() {
    return !isEmpty() && isClippingOperator(cigarElements[numCigarElements()-1].getOperator());
}

bool Cigar::isClipped() {
    return isLeftClipped() || isRightClipped();
}

std::vector<std::logic_error> Cigar::isValid(std::string readName, long recordNumber) {
    if(isEmpty()) {
        return {};
    }
    std::vector<std::logic_error> ret;
    bool seenRealOperator = false;
    for(int i = 0; i < cigarElements.size(); ++i) {
        CigarElement element = cigarElements[i];
        if(element.getLength() == 0) {
            ret.emplace_back(std::logic_error("CIGAR element with zero length"));
        }
        CigarOperator op = element.getOperator();

        if(isClippingOperator(op)) {
            if(op == H) {
                if(i != 0 && i != cigarElements.size() - 1) {
                    ret.emplace_back(std::logic_error("Hard clipping operator not at start or end of CIGAR"));
                }
            }
            else {
                if(op != S)
                    throw std::invalid_argument("Should never happen");
                if(i == 0 || i == cigarElements.size() - 1) {

                }
                else if (i == 1) {
                    if (cigarElements.size() == 3 && cigarElements[2].getOperator() == H) {

                    } else if (cigarElements[0].getOperator() != H) {
                        ret.emplace_back(std::logic_error(
                                "Soft clipping CIGAR operator can only be inside of hard clipping operator"));
                    }
                }
                else if (i == cigarElements.size() - 2) {
                    if(cigarElements[cigarElements.size() - 1].getOperator() != H) {
                        ret.emplace_back(std::logic_error("Soft clipping CIGAR operator can only be inside of hard clipping operator"));
                    }
                }
                else {
                    ret.emplace_back(std::logic_error("Soft clipping CIGAR operator can at start or end of read, or be inside of hard clipping operator"));
                }
            }
        }
        else if (isRealOperator(op)) {
            seenRealOperator = true;
            if(isInDelOperator(op)) {
                for(int j = i+1; j < cigarElements.size(); ++j) {
                    CigarOperator nextOperator = cigarElements[j].getOperator();
                    if ((isRealOperator(nextOperator) && !isInDelOperator(nextOperator)) || isPaddingOperator(nextOperator)) {
                        break;
                    }
                    if (isInDelOperator(nextOperator) && op == nextOperator) {
                        ret.emplace_back(std::logic_error("No M or N operator between"));
                    }
                }
            }
        }
        else if (isPaddingOperator(op)) {
            if(i == 0) {

            }
            else if (i == cigarElements.size() - 1) {
                ret.emplace_back(std::logic_error("Padding operator not valid at end of CIGAR"));
            }
            else if (!isRealOperator(cigarElements[i-1].getOperator()) ||
                    !isRealOperator(cigarElements[i+1].getOperator())) {
                ret.emplace_back(std::logic_error("Padding operator not between real operators in CIGAR"));
            }
        }
    }
    if(!seenRealOperator) {
        ret.emplace_back(std::logic_error("No real operator (M|I|D|N) in CIGA"));
    }
    return ret;
}

bool Cigar::isValid(bam1_t *read)
{
    if(read->core.n_cigar == 0)
        return false;

    bool seenRealOperator = false;
    uint32_t * cigarArray = bam_get_cigar(read);
    for(int i=0; i<read->core.n_cigar; i++){
        uint32_t cigar = cigarArray[i];
        if(bam_cigar_oplen(cigar) == 0)
            return false;

        if(bam_cigar_op(cigar) == BAM_CHARD_CLIP)
        {
            if(i != 0 && i != read->core.n_cigar-1)
                return false;
        }
        else if(bam_cigar_op(cigar) == BAM_CSOFT_CLIP) {
            if (i == 0 || i == read->core.n_cigar - 1) {
                // Soft clip at either end is fine
            } else if (i == 1) {
                if(read->core.n_cigar == 3 && bam_cigar_op(cigarArray[2]) == BAM_CHARD_CLIP){

                } else if(bam_cigar_op(cigarArray[0]) != BAM_CHARD_CLIP){
                    return false;
                }
            } else if(i == read->core.n_cigar - 2) {
                if(bam_cigar_op(cigarArray[read->core.n_cigar-1]) != BAM_CHARD_CLIP){
                    return false;
                }
            } else {
                return false;
            }
        } else if(isRealOperator(bam_cigar_op(cigar))) {
            // Must be at least one real operator (MIDN)
            seenRealOperator = true;
            // There should be an M or P operator between any pair of IDN operators
            if (isInDelOperator(bam_cigar_op(cigar))) {
                for (int j = i+1; j < read->core.n_cigar; ++j) {
                    uint32_t nextOperator = bam_cigar_op( cigarArray[j]);
                    // Allow
                    if ((isRealOperator(nextOperator) && !isInDelOperator(nextOperator)) || isPaddingOperator(nextOperator)) {
                        break;
                    }
                    if (isInDelOperator(nextOperator) && bam_cigar_op(cigar) == nextOperator) {
                        return false;
                    }
                }
            }
        } else if(isPaddingOperator(bam_cigar_op(cigar))){
            if (i == 0) {
                /*
                 * Removed restriction that padding not be the first operator because if a read starts in the middle of a pad
                 * in a padded reference, it is necessary to precede the read with padding so that alignment start refers to a
                 * position on the unpadded reference.
                */
            } else if (i == read->core.n_cigar - 1) {
                return false;
            } else if (!isRealOperator(bam_cigar_op(cigarArray[i-1])) || !isRealOperator(bam_cigar_op(cigarArray[i+1]))) {
                return false;
            }
        }
    }

    if(!seenRealOperator){
        return false;
    }
    return true;
}

CigarElement Cigar::getFirstCigarElement() {
    return isEmpty() ? CigarElement(0, M) : cigarElements.at(0);
}

CigarElement Cigar::getLastCigarElement() {
    return isEmpty() ? CigarElement(0, M) : cigarElements.at(cigarElements.size() - 1);
}
