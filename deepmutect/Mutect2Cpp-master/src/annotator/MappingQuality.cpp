//
// Created by lhh on 6/15/22.
//

#include "MappingQuality.h"

int MappingQuality::aggregate(std::vector<int> &values) {
    return values.empty() ? VALUE_FOR_NO_READS : MathUtils::median(values);
}

std::string MappingQuality::getVcfKey() {
    return VCFConstants::MEDIAN_MAPPING_QUALITY_KEY;
}

std::string MappingQuality::getDescription() {
    return "median mapping quality";
}

std::optional<int> MappingQuality::getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc) {
    assert(read != nullptr);
    return {read->getMappingQuality()};
}