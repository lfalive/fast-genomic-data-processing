//
// Created by 梦想家xixi on 2021/12/21.
//

#include "ClippingOp.h"
#include <memory>
#include <stack>
#include <utility>
#include "ReadUtils.h"

ClippingOp::ClippingOp(int start, int stop) : start(start), stop(stop){
}

std::shared_ptr<SAMRecord> ClippingOp::applyHardClipBases(std::shared_ptr<SAMRecord> read, int start, int stop) {
    std::shared_ptr<Cigar> cigar = read->getCigar();
    std::shared_ptr<Cigar> tmp = nullptr;
    std::shared_ptr<ClippingOp::CigarShift> cigarShift;
    if(read->isUnmapped()) {
        tmp = std::make_shared<Cigar>();
        cigarShift = std::make_shared<ClippingOp::CigarShift>(tmp, 0, 0);
    } else {
        cigarShift = hardClipCigar(cigar, start, stop);
    }
    int newLength = read->getLength() - (stop - start + 1) - cigarShift->shiftFromStart - cigarShift->shiftFromEnd;
    if(newLength == 0) {
        return ReadUtils::emptyRead(read);
    }
    std::shared_ptr<uint8_t[]> newBases(new uint8_t[newLength]{0});
    std::shared_ptr<uint8_t[]> newQuals(new uint8_t[newLength]{0});
    int copyStart = (start == 0) ? stop + 1 + cigarShift->shiftFromStart : cigarShift->shiftFromStart;
    memcpy(newBases.get(), read->getBasesNoCopy().get()+copyStart, newLength);
    memcpy(newQuals.get(), read->getBaseQualitiesNoCopy().get()+copyStart, newLength);
    std::shared_ptr<SAMRecord> hardClippedRead{new SAMRecord(*read)};
    hardClippedRead->setBaseQualities(newQuals, newLength);
    hardClippedRead->setBases(newBases, newLength);
    hardClippedRead->setCigar(cigarShift->cigar);
    if(start == 0 && !read->isUnmapped()) {
        hardClippedRead->setPosition(read->getContig(), read->getStart() + calculateAlignmentStartShift(cigar, cigarShift->cigar));
    }
    if(ReadUtils::hasBaseIndelQualities(read)) {
        std::shared_ptr<uint8_t[]> newBaseInsertionQuals(new uint8_t[newLength]);
        std::shared_ptr<uint8_t[]> newBaseDeletionQuals(new uint8_t[newLength]);
        int length;
        std::shared_ptr<uint8_t[]> new_base = ReadUtils::getBaseInsertionQualities(read, length);

        memcpy(newBaseInsertionQuals.get(), new_base.get()+copyStart, newLength);
        new_base = ReadUtils::getBaseDeletionQualities(read, length);
        memcpy(newBaseDeletionQuals.get(), new_base.get()+copyStart, newLength);
        ReadUtils::setDeletionBaseQualities(hardClippedRead, newBaseDeletionQuals, newLength);
        ReadUtils::setInsertionBaseQualities(hardClippedRead, newBaseInsertionQuals, newLength);
    }
    return hardClippedRead;
}

std::shared_ptr<ClippingOp::CigarShift> ClippingOp::hardClipCigar(std::shared_ptr<Cigar> cigar, int start, int stop) {
    std::shared_ptr<Cigar> newCigar(new Cigar());
    int index = 0;
    int totalHardClipCount = stop - start + 1;
    int alignmentShift = 0;

    if(start == 0) {
        std::vector<CigarElement> cigarElementIterator = cigar->getCigarElements();
        int i = 0;
        CigarElement cigarElement = cigarElementIterator[0];
        i++;
        while(cigarElement.getOperator() == H) {
            totalHardClipCount += cigarElement.getLength();
            if(i < cigarElementIterator.size()) {
                cigarElement = cigarElementIterator[i];
                i++;
            } else {
                throw std::invalid_argument("Read is entirely hard-clipped, shouldn't be trying to clip it's cigar string");
            }
        }
        while(index <= stop) {
            int shift = 0;
            if(CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                shift = cigarElement.getLength();
            }
            if(index + shift == stop + 1) {
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());
                newCigar->add(CigarElement(totalHardClipCount + alignmentShift, H));
            } else if (index + shift > stop + 1) {
                int elementLengthAfterChopping = cigarElement.getLength() - (stop - index + 1);
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, stop - index + 1);
                newCigar->add(CigarElement(totalHardClipCount + alignmentShift, H));
                newCigar->add(CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
            }
            index += shift;
            alignmentShift += calculateHardClippingAlignmentShift(cigarElement, shift);

            if(index <= stop && i < cigarElementIterator.size()) {
                cigarElement = cigarElementIterator[i];
                i++;
            } else {
                break;
            }
        }

        while(i < cigarElementIterator.size()) {
            cigarElement = cigarElementIterator[i];
            newCigar->add(CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
            i++;
        }
    } else {
        std::vector<CigarElement> cigarElementIterator = cigar->getCigarElements();
        int i = 0;
        CigarElement cigarElement = cigarElementIterator[i];
        i++;

        while(index < start) {
            int shift = 0;
            if(CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
                shift = cigarElement.getLength();
            }

            if(index + shift < start) {
                newCigar->add(CigarElement(cigarElement.getLength(), cigarElement.getOperator()));
            } else {
                int elementLengthAfterChopping = start - index;
                alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength() - (start - index));

                if(cigarElement.getOperator() == H) {
                    totalHardClipCount += elementLengthAfterChopping;
                } else {
                    newCigar->add(CigarElement(elementLengthAfterChopping, cigarElement.getOperator()));
                }
            }
            index += shift;
            if(index < start && i < cigarElementIterator.size()) {
                cigarElement = cigarElementIterator[i];
                i++;
            } else {
                break;
            }
        }
        while(i < cigarElementIterator.size()) {
            cigarElement = cigarElementIterator[i];
            i++;
            alignmentShift += calculateHardClippingAlignmentShift(cigarElement, cigarElement.getLength());

            if(cigarElement.getOperator() == H) {
                totalHardClipCount += cigarElement.getLength();
            }
        }
        newCigar->add(CigarElement(totalHardClipCount + alignmentShift, H));
    }
    return cleanHardClippedCigar(newCigar);
}

int ClippingOp::calculateHardClippingAlignmentShift(CigarElement &cigarElement, int clippedLength) {
    if(cigarElement.getOperator() == I) {
        return -clippedLength;
    } else if (cigarElement.getOperator() == D || cigarElement.getOperator() == N) {
        return cigarElement.getLength();
    }
    return 0;
}

std::shared_ptr<ClippingOp::CigarShift> ClippingOp::cleanHardClippedCigar(std::shared_ptr<Cigar> cigar) {
    std::shared_ptr<Cigar> cleanCigar(new Cigar());
    int shiftFromStart = 0;
    int shiftFromEnd = 0;
    std::stack<CigarElement> cigarStack;
    std::stack<CigarElement> inverseCigarStack;
    for(CigarElement cigarElement : cigar->getCigarElements()) {
        cigarStack.push(cigarElement);
    }
    for(int i = 1; i <= 2; i++) {
        int shift = 0;
        int totalHardClip = 0;
        bool readHasStarted = false;
        bool addedHardClips = false;

        while(!cigarStack.empty()) {
            CigarElement cigarElement = cigarStack.top();
            cigarStack.pop();
            if(!readHasStarted &&
            cigarElement.getOperator() != D &&
            cigarElement.getOperator() != N &&
            cigarElement.getOperator() != H) {
                readHasStarted = true;
            } else if (!readHasStarted && cigarElement.getOperator() == H) {
                totalHardClip += cigarElement.getLength();
            } else if (!readHasStarted && cigarElement.getOperator() == D) {
                totalHardClip += cigarElement.getLength();
            } else if (!readHasStarted && cigarElement.getOperator() == N) {
                totalHardClip += cigarElement.getLength();
            }

            if(readHasStarted) {
                if(i == 1) {
                    if(!addedHardClips) {
                        if(totalHardClip > 0) {
                            inverseCigarStack.push(CigarElement(totalHardClip, H));
                        }
                        addedHardClips = true;
                    }
                    inverseCigarStack.push(cigarElement);
                } else {
                    if(!addedHardClips) {
                        if(totalHardClip > 0) {
                            cleanCigar->add(CigarElement(totalHardClip, H));
                        }
                        addedHardClips = true;
                    }
                    cleanCigar->add(cigarElement);
                }
            }
        }
        if(i == 1) {
            shiftFromEnd = shift;
            cigarStack = inverseCigarStack;
        } else {
            shiftFromStart = shift;
        }
    }
    return std::make_shared<ClippingOp::CigarShift>(cleanCigar, shiftFromStart, shiftFromEnd);
}

int ClippingOp::calculateAlignmentStartShift(std::shared_ptr<Cigar> oldCigar, std::shared_ptr<Cigar> newCigar) {
    int newShift = calcHardSoftOffset(newCigar);
    int oldShift = calcHardSoftOffset(oldCigar);
    return newShift - oldShift;
}

int ClippingOp::calcHardSoftOffset(std::shared_ptr<Cigar> cigar) {
    std::vector<CigarElement> elements = cigar->getCigarElements();
    int size = 0;
    int i = 0;
    while ( i < elements.size() && elements[i].getOperator() == H ) {
        size += elements[i].getLength();
        i++;
    }
    while ( i < elements.size() && elements[i].getOperator() == S ) {
        size += elements[i].getLength();
        i++;
    }
    return size;
}

std::shared_ptr<SAMRecord> ClippingOp::apply(ClippingRepresentation algorithm, std::shared_ptr<SAMRecord> originalRead, bool runAsserts) {
    switch(algorithm){
        case HARDCLIP_BASES: {
            return applyHardClipBases(originalRead, start, stop);
        }
        case REVERT_SOFTCLIPPED_BASES: {
            return applyRevertSoftClippedBases(originalRead);
        }
        default: {
            throw std::invalid_argument("Unexpected Clipping operator type");
        }
    }
}

std::shared_ptr<SAMRecord> ClippingOp::applyRevertSoftClippedBases(const std::shared_ptr<SAMRecord>& read) {
    std::shared_ptr<SAMRecord> unclipped(new SAMRecord(*read));
    std::shared_ptr<Cigar> unclippedCigar(new Cigar());
    int matchesCount = 0;
    for(CigarElement element : read->getCigarElements()) {
        if(element.getOperator() == S || element.getOperator() == M) {
            matchesCount += element.getLength();
        }else if (matchesCount > 0) {
            unclippedCigar->add(CigarElement(matchesCount, M));
            matchesCount = 0;
            unclippedCigar->add(element);
        } else {
            unclippedCigar->add(element);
        }
    }
    if (matchesCount > 0) {
        unclippedCigar->add(CigarElement(matchesCount, M));
    }
    unclipped->setCigar(unclippedCigar);
    int newStart = read->getStart() + calculateAlignmentStartShift(read->getCigar(), unclippedCigar);
    if(newStart < 0) {
        unclipped->setPosition(unclipped->getContig(), 0);
        unclipped = applyHardClipBases(unclipped, 0, (-newStart) - 1);

        if(! unclipped->isUnmapped()) {
            unclipped->setPosition(unclipped->getContig(), 0);
        }
        return unclipped;
    }
    else {
        unclipped->setPosition(unclipped->getContig(), newStart);
        return unclipped;
    }
}

ClippingOp::CigarShift::CigarShift(std::shared_ptr<Cigar> cigar, int shiftFromStart, int shiftFromEnd) : cigar(std::move(cigar)), shiftFromStart(shiftFromStart), shiftFromEnd(shiftFromEnd){
}
