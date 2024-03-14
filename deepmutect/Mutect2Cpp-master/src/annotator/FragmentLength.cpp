//
// Created by lhh on 6/14/22.
//

#include "FragmentLength.h"
#include <optional>

int FragmentLength::aggregate(std::vector<int> &values) {
    return values.empty() ? 0 : MathUtils::median(values);
}

bool FragmentLength::includeRefAllele() {
    return true;
}

std::string FragmentLength::getVcfKey() {
    return VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY;
}

std::string FragmentLength::getDescription() {
    return  "median fragment length";
}

std::optional<int> FragmentLength::getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc) {
    assert(read != nullptr);

    //abs because fragment lengths are negative if mate comes first
    return {abs(read->getFragmentLength())};
}