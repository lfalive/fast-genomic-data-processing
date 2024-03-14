//
// Created by 梦想家xixi on 2021/10/12.
//

#include "Locatable.h"


int Locatable::getLengthOnReference() {
    return this->getEnd() - this->getStart() + 1;
}

bool Locatable::contigsMatch(const std::shared_ptr<Locatable> &other) {
    return (!this->getContig().empty()) && other != nullptr && this->getContig() == other->getContig();
}

bool Locatable::contains(const std::shared_ptr<Locatable> &other) {
    return this->contigsMatch(other) && Mutect2Utils::encloses(this->getStart(), this->getEnd(), other->getStart(), other->getEnd());
}

bool Locatable::withinDistanceOf(const std::shared_ptr<Locatable> &other, int distance) {
    return this->contigsMatch(other) && Mutect2Utils::overlaps(this->getStart(), this->getEnd(), other->getStart() - distance, other->getEnd() + distance);
}

bool Locatable::overlaps(const std::shared_ptr<Locatable> & other) {
    return this->withinDistanceOf(other, 0);
}