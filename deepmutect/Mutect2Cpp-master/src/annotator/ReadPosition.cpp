//
// Created by lhh on 6/15/22.
//

#include "ReadPosition.h"
#include "ReadPosRankSumTest.h"

int ReadPosition::aggregate(std::vector<int> &values) {
    return values.empty() ? VALUE_FOR_NO_READS : MathUtils::median(values);
}

std::string ReadPosition::getVcfKey() {
    return VCFConstants::MEDIAN_READ_POSITON_KEY;
}

std::string ReadPosition::getDescription() {
    return "median distance from end of read";
}

std::optional<int> ReadPosition::getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc) {
    if (vc->getStart() < read->getStart() || read->getEnd() < vc->getStart()) {
        return std::nullopt;
    }

    auto valueAsDouble = ReadPosRankSumTest::getReadPosition(read, vc->getStart());
    return valueAsDouble.has_value() ? std::optional<int>((int)valueAsDouble.value()) : std::nullopt;
}