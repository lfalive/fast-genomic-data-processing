//
// Created by lhh on 6/2/22.
//

#include "BaseQuality.h"
#include "BaseQualityRankSumTest.h"

int BaseQuality::aggregate(std::vector<int> &values) {
    return values.empty() ? 0 : MathUtils::median(values);
}

std::string BaseQuality::getVcfKey() {
    return VCFConstants::MEDIAN_BASE_QUALITY_KEY;
}

std::string BaseQuality::getDescription() {
    return "median base quality";
}

bool BaseQuality::includeRefAllele() {
    return true;
}

std::optional<int> BaseQuality::getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc) {
    return getBaseQuality(read, vc);
}

std::optional<int> BaseQuality::getBaseQuality(const std::shared_ptr<SAMRecord>& read, const shared_ptr<VariantContext>& vc) {
    if (vc->getStart() < read->getStart() || read->getEnd() < vc->getStart()) {
        return std::nullopt;
    }

    auto result = BaseQualityRankSumTest::getReadBaseQuality(read, vc->getStart());
    return result.has_value() ? std::optional<int>(round(result.value())) : std::nullopt;
}