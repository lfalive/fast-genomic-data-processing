//
// Created by cluster on 22-11-8.
//

#include "SequencingError.h"
#include <stdexcept>

double SequencingError::logLikelihood(const Datum& datum) {
    return 0;
}

double SequencingError::logLikelihood(int totalCount, int altCount) {
    throw std::invalid_argument("This method should never be called on the sequencing error cluster.");
}

void SequencingError::learn(std::vector<Datum> data) {

}

std::string SequencingError::toString() {
    return ("sequencing error");
}

SequencingError::SequencingError() = default;

SequencingError::~SequencingError() = default;
