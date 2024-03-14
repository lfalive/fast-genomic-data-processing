//
// Created by cluster on 22-11-8.
//

#include "BetaDistributionShape.h"
#include <stdexcept>

BetaDistributionShape BetaDistributionShape::FLAT_BETA = BetaDistributionShape(1, 1);

BetaDistributionShape::BetaDistributionShape(double alpha, double beta) : alpha(alpha), beta(beta) {
    if(alpha < 0 || beta < 0) {
        std::string param = ("alpha and beta must be greater than 0 but got");
        throw std::invalid_argument(param);
    }
}

double BetaDistributionShape::getAlpha() const{
    return alpha;
}

double BetaDistributionShape::getBeta() const{
    return beta;
}

void BetaDistributionShape::initial() {
    BetaDistributionShape(1, 1);
}
