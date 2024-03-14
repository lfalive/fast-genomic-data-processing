//
// Created by lhh on 5/24/22.
//

#include "parallel_hashmap/phmap.h"
#include "VariantContextUtils.h"

int VariantContextUtils::getSize(std::shared_ptr<VariantContext> &vc)
{
    return vc->getEnd() - vc->getStart() + 1;
}

void VariantContextUtils::calculateChromosomeCounts(std::shared_ptr<VariantContext> vc, std::map<std::string,std::string> &attributes,
                                                    bool removeStaleValues) {
    phmap::flat_hash_set<std::string> founderIds;
    int AN = vc->getCalledChrCount();

    // if everyone is a no-call, remove the old attributes if requested
    if ( AN == 0 && removeStaleValues ) {
        if ( attributes.find(VCFConstants::ALLELE_COUNT_KEY) != attributes.end())
            attributes.erase(VCFConstants::ALLELE_COUNT_KEY);
        if ( attributes.find(VCFConstants::ALLELE_FREQUENCY_KEY) != attributes.end())
            attributes.erase(VCFConstants::ALLELE_FREQUENCY_KEY);
        if ( attributes.find(VCFConstants::ALLELE_NUMBER_KEY)  != attributes.end() )
            attributes.erase(VCFConstants::ALLELE_NUMBER_KEY);
        return;
    }

    if(vc->hasGenotypes())
    {
        //attributes.insert(VCFConstants::ALLELE_NUMBER_KEY, AN); // TODO: how to represent attribute?

    }

}