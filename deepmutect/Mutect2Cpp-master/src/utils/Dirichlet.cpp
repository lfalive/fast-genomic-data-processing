//
// Created by lhh on 5/28/22.
//

#include <cassert>
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
#include "Dirichlet.h"
#include "Mutect2Utils.h"
#include "MathUtils.h"

Dirichlet::Dirichlet(std::shared_ptr<std::vector<double>> alpha) {
    assert(alpha != nullptr);
    Mutect2Utils::validateArg(alpha->size() >= 1, "Dirichlet parameters must have at least one element");
    this->alpha.reserve(alpha->size());
    for(double a : *alpha)
    {
        Mutect2Utils::validateArg(a >= 0 && !std::isinf(a), "Dirichlet parameters may not be negative or infinite");
        this->alpha.push_back(a);
    }
}

std::shared_ptr<std::vector<double>> Dirichlet::effectiveLogMultinomialWeights() {
    double digammaOfSum = boost::math::digamma(MathUtils::sum(alpha));
    auto result = make_shared<vector<double>>(alpha.size());
    assert(result->size() == alpha.size());
    for(int m=0; m<result->size(); m++)
    {
        result->operator[](m) = boost::math::digamma(alpha[m]) - digammaOfSum;
    }
    return result;
}