//
// Created by lhh on 6/13/22.
//

#include "Coverage.h"

std::shared_ptr<std::map<std::string, AttributeValue>>
Coverage::annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc,
                   AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);
    if (likelihoods == nullptr || likelihoods->evidenceCount() == 0) {
        return nullptr;
    }

    int depth = likelihoods->evidenceCount();
    auto result =  make_shared<std::map<std::string, AttributeValue>>();
    result->emplace(getKeyNames().at(0), AttributeValue(depth));
    return result;
}

std::vector<std::string> Coverage::getKeyNames() {
    return {VCFConstants::DEPTH_KEY};
}