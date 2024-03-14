//
// Created by lhh on 6/2/22.
//

#include "PerAlleleAnnotation.h"

std::shared_ptr<std::map<std::string, AttributeValue>> PerAlleleAnnotation::annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc,
                                   AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);
    if(!likelihoods)
        return nullptr;

    std::unordered_map<Allele*, vector<int>> values;
    for(auto& allele: likelihoods->getAlleles())
    {
        values.emplace(allele.get(), vector<int>());
    }

    auto bestAlleles = likelihoods->bestAllelesBreakingTies();
    for(const auto& ba: *bestAlleles)
    {
        if(ba->isInformative() && isUsableRead(ba->evidence))
        {
            auto value = getValueForRead(ba->evidence, vc);
            if(value.has_value())
                values.at(ba->allele.get()).emplace_back(value.value());
        }
    }

    vector<int> statistics;
    for(auto& allele : vc->getAlleles())
    {
        if(includeAllele(allele))
        {
            statistics.emplace_back(aggregate(values.at(allele.get())));
        }
    }

    auto result = make_shared<std::map<std::string, AttributeValue>>();
    result->emplace(getVcfKey(), AttributeValue(statistics));
    return result;
}

bool PerAlleleAnnotation::isUsableRead(std::shared_ptr<SAMRecord> read) {
    return read->getMappingQuality() != 0 && read->getMappingQuality() != QualityUtils::MAPPING_QUALITY_UNAVALIABLE;
}

bool PerAlleleAnnotation::includeAllele(shared_ptr<Allele> allele) {
    return allele->getIsNonReference() || includeRefAllele();
}

std::vector<std::string> PerAlleleAnnotation::getKeyNames() {
    return {getVcfKey()};
}
