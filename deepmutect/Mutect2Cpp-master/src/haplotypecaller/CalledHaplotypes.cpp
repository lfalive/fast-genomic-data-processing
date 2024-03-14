//
// Created by lhh on 5/20/22.
//

#include <cassert>
#include <memory>

#include "CalledHaplotypes.h"


CalledHaplotypes::CalledHaplotypes(std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> calls, std::shared_ptr<phmap::flat_hash_set<std::shared_ptr<Haplotype>>> calledHaplotypes) : calls(calls), calledHaplotypes(calledHaplotypes)
{
    assert(calls->empty() == calledHaplotypes->empty());

}

std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> CalledHaplotypes::getCalls() {
    return calls;
}

std::shared_ptr<phmap::flat_hash_set<std::shared_ptr<Haplotype>>> CalledHaplotypes::getCalledHaplotypes() {
    return calledHaplotypes;
}