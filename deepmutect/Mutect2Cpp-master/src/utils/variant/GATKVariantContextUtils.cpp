//
// Created by lhh on 4/26/22.
//

#include <cassert>
#include <utils/Utils.h>
#include "parallel_hashmap/phmap.h"
#include <algorithm>
#include "parallel_hashmap/phmap.h"
#include "GATKVariantContextUtils.h"
#include "VariantContextUtils.h"
#include "variantcontext/builder/GenotypeBuilder.h"
#include "variantcontext/builder/VariantContextBuilder.h"
#include "Mutect2Utils.h"

AlleleMapper::AlleleMapper(std::shared_ptr<VariantContext> vc): vc(vc) {

}

AlleleMapper::AlleleMapper(std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>> map): map(map){

}

bool AlleleMapper::needsRemapping() {
    return this->map != nullptr;
}

std::shared_ptr<Allele> AlleleMapper::remap(std::shared_ptr<Allele> &a)
{
    return map != nullptr && (map->find(a) != map->end()) ? map->at(a) : a;
}

std::shared_ptr<std::vector<std::shared_ptr<Allele>>> AlleleMapper::remap(std::vector<std::shared_ptr<Allele>> &as)
{
    std::shared_ptr<std::vector<std::shared_ptr<Allele>>> newAs(new std::vector<std::shared_ptr<Allele>>);
    for(auto& a : as)
    {
        newAs->emplace_back(remap(a));
    }
    return newAs;
}

std::shared_ptr<std::vector<std::shared_ptr<Allele>>> AlleleMapper::values()
{
    std::shared_ptr<std::vector<std::shared_ptr<Allele>>> alleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
    if(map == nullptr || map->empty())
    {
        *alleles = vc->getAlleles();
    } else {
        alleles->reserve(map->size());
        for(auto & iter : *map)
        {
            alleles->emplace_back(iter.second);
        }
    }
    return alleles;
}

int GATKVariantContextUtils::findNumberOfRepetitions(uint8_t *repeatUnitFull, int repeatUnitFullLength, int offsetInRepeatUnitFull,
                                                     int repeatUnitLength, uint8_t *testStringFull, int testStringFullLength,
                                                     int offsetInTestStringFull, int testStringLength,
                                                     bool leadingRepeats) {
    if (testStringLength == 0){
        return 0;
    }

    assert(repeatUnitLength >= 0 && repeatUnitLength <= repeatUnitFullLength);
    assert(offsetInRepeatUnitFull >= 0 && offsetInRepeatUnitFull < repeatUnitFullLength);
    assert(offsetInTestStringFull >= 0 && offsetInTestStringFull < testStringFullLength);
    assert(testStringLength >= 0 && testStringLength <= testStringFullLength);

    int lengthDifference = testStringLength - repeatUnitLength;

    if(leadingRepeats)
    {
        int numRepeats = 0;
        // look forward on the test string
        for (int start = 0; start <= lengthDifference; start += repeatUnitLength) {
            if(Utils::equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                numRepeats++;
            } else {
                return numRepeats;
            }
        }
        return numRepeats;
    } else {
        // look backward. For example, if repeatUnit = AT and testString = GATAT, number of repeat units is still 2
        int numRepeats = 0;
        // look backward on the test string
        for (int start = lengthDifference; start >= 0; start -= repeatUnitLength) {
            if (Utils::equalRange(testStringFull, start + offsetInTestStringFull, repeatUnitFull, offsetInRepeatUnitFull, repeatUnitLength)) {
                numRepeats++;
            } else {
                return numRepeats;
            }
        }
        return numRepeats;
    }

}

int GATKVariantContextUtils::findNumberOfRepetitions(uint8_t* repeatUnit, int repeatUnitLength, uint8_t* testString, int testStringLength, bool leadingRepeats) {
    if(testStringLength == 0)
        return 0;
    return findNumberOfRepetitions(repeatUnit, repeatUnitLength, 0, repeatUnitLength, testString, testStringLength, 0, testStringLength, leadingRepeats);
}

std::shared_ptr<VariantContext> GATKVariantContextUtils::simpleMerge(std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> unsortedVCs,
                                                                     std::vector<std::string>& priorityListOfVCs,
                                                                     FilteredRecordMergeType filteredRecordMergeType,
                                                                     GenotypeMergeType genotypeMergeOptions,
                                                                     bool filteredAreUncalled) {
    if(unsortedVCs == nullptr || unsortedVCs->empty())
        return nullptr;

    int originalNumOfVCs = priorityListOfVCs.size();
    auto preFilteredVCs = sortVariantContextsByPriority(unsortedVCs, priorityListOfVCs, genotypeMergeOptions);
    // Make sure all variant contexts are padded with reference base in case of indels if necessary
    std::vector<std::shared_ptr<VariantContext>> VCs;
    for(auto vc : *preFilteredVCs)
    {
        if(!filteredAreUncalled || vc->isNotFiltered())
            VCs.emplace_back(vc);
    }

    if ( VCs.empty() ) // everything is filtered out and we're filteredAreUncalled
        return nullptr;

    // establish the baseline info from the first VC
    std::shared_ptr<VariantContext> first = VCs[0];
    std::string name = first->getSource();
    auto refAllele = determineReferenceAllele(VCs);

    std::shared_ptr<std::vector<std::shared_ptr<Allele>>> alleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();   //---use a vector to function as a set
    std::set<std::string> filters;
    std::map<std::string, std::string> attributes;
    std::set<std::string> inconsistentAttributes;
    std::set<std::string> variantSources;  // contains the set of sources we found in our set of VCs that are variant
    std::set<std::string> rsIDs;

    std::shared_ptr<VariantContext> longestVC = first;
    int depth = 0;
    double log10PError = CommonInfo::NO_LOG10_PERROR;
    bool anyVCHadFiltersApplied = false;
    auto genotypes = std::make_shared<GenoTypesContext>();

    // counting the number of filtered and variant VCs
    int nFiltered = 0;

    // cycle through and add info from the other VCs, making sure the loc/reference matches
    for(auto& vc : VCs)
    {
        assert(longestVC->getStart() == vc->getStart());

        if(VariantContextUtils::getSize(vc) > VariantContextUtils::getSize(longestVC))
            longestVC = vc;

        nFiltered += vc->isFiltered() ? 1 : 0;
        if(vc->isVariant())
            variantSources.insert(vc->getSource());

        auto alleleMapping = resolveIncompatibleAlleles(refAllele, vc, *alleles);

        auto temp = alleleMapping->values();
        for(auto value : *temp)
        {
            bool present = false;
            for(auto a : *alleles)
            {
                if(*a == *value)
                {
                    present = true;
                    break;
                }
            }
            if(!present)
                alleles->emplace_back(value);
        }

        mergeGenotypes(*genotypes, vc, alleleMapping, genotypeMergeOptions == UNIQUIFY);

        // We always take the QUAL of the first VC with a non-MISSING qual for the combined value
        if(log10PError == CommonInfo::NO_LOG10_PERROR)
            log10PError = vc->getLog10PError();

        for(auto& filter : vc->getFilter())
            filters.insert(filter);
        anyVCHadFiltersApplied |= vc->filtersWereApplied();

        //
        // add attributes
        //
        // special case DP (add it up) and ID (just preserve it)
        //
        if(vc->hasAttribute(VCFConstants::DEPTH_KEY))
            depth += vc->getAttributeAsInt(VCFConstants::DEPTH_KEY, 0);

        if(vc->hasID())
            rsIDs.insert(vc->getID());

        for(auto& iter : vc->getAttributes())
        {
            std::string key = iter.first;
            AttributeValue value = iter.second;

            // only output annotations that have the same value in every input VC
            // if we don't like the key already, don't go anywhere
            // TODO: finish this method
           /* if(inconsistentAttributes.find(key) == inconsistentAttributes.end())
            {
                bool alreadyFound = attributes.find(key) != attributes.end();
                void * boundValue = attributes.at(key);
                bool boundIsMissingValue = alreadyFound && boundValue == nullptr;    //---in GATK, this line is: boolean boundIsMissingValue = alreadyFound && boundValue.equals(VCFConstants.MISSING_VALUE_v4);

                if(alreadyFound && !(*boundValue == *value) && !boundIsMissingValue)
                {

                }
            }*/
        }
    }

    // if we have more alternate alleles in the merged VC than in one or more of the
    // original VCs, we need to strip out the GL/PLs (because they are no longer accurate), as well as allele-dependent attributes like AC,AF, and AD
    for(auto& vc : VCs)
    {
        if(vc->getAlleles().size() == 1)
            continue;


        if(hasPLIncompatibleAlleles(*alleles, vc->getAlleles()))
        {
            if(!genotypes->isEmpty())
            {
                throw std::exception();
            }
            auto temp = stripPLsAndAD(genotypes);
            genotypes = temp;

            // this will remove stale AC,AF attributed from vc
            VariantContextUtils::calculateChromosomeCounts(vc, attributes, true);
            break;
        }
    }

    // if at least one record was unfiltered and we want a union, clear all of the filters
    if((filteredRecordMergeType == KEEP_IF_ANY_UNFILTERED && nFiltered != VCs.size()) || filteredRecordMergeType == KEEP_UNCONDITIONAL)
        filters.clear();

    if(depth > 0)
        attributes.insert({VCFConstants::DEPTH_KEY, std::to_string(depth)});

    std::string ID = rsIDs.empty() ? VCFConstants::EMPTY_ID_FIELD : Utils::join(",", rsIDs);

    auto builder = std::shared_ptr<VariantContextBuilder>((new VariantContextBuilder())->setSource(name)->setId(ID));
    builder->setLoc(longestVC->getContig(), longestVC->getStart(), longestVC->getEnd());
    builder->setAlleles(alleles);
    builder->setGenotypes(genotypes);
    builder->setLog10PError(log10PError);
    if(anyVCHadFiltersApplied){
        builder->setFilters(new std::set<std::string>);
    }

    // Trim the padded bases of all alleles if necessary
    auto merged = builder->make();
    return merged;
}

std::shared_ptr<GenoTypesContext> GATKVariantContextUtils::stripPLsAndAD(std::shared_ptr<GenoTypesContext> genotypes)
{
    auto newGs = std::make_shared<GenoTypesContext>(genotypes->getSize());
    for(int i=0; i<genotypes->getSize(); i++)
    {
        newGs->add(removePLsAndAD(genotypes->get(i)));
    }
    return newGs;
}

std::shared_ptr<Genotype> GATKVariantContextUtils::removePLsAndAD(std::shared_ptr<Genotype> g) {
    return (g->hasLikelihoods() || g->hasAD()) ? GenotypeBuilder(g).make() : g;
}

bool GATKVariantContextUtils::hasPLIncompatibleAlleles(std::vector<std::shared_ptr<Allele>> &alleleSet1,
                                                       std::vector<std::shared_ptr<Allele>> &alleleSet2) {
    auto it1 = alleleSet1.begin();
    auto it2 = alleleSet2.begin();

    while(it1 != alleleSet1.end() && it2 != alleleSet2.end())
    {
        if(!(*it1 == *it2 || it1->get() == it2->get() || *(it1->get()) == *(it2->get())))
        {
            return true;
        }
        it1++;
        it2++;
    }

    // by this point, at least one of the iterators is empty.  All of the elements
    // we've compared are equal up until this point.  But it's possible that the
    // sets aren't the same size, which is indicated by the test below.  If they
    // are of the same size, though, the sets are compatible
    return it1 != alleleSet1.end() || it2 != alleleSet2.end();
}

class CompareByPriority{
    phmap::flat_hash_map<std::string, int>& ComparatorMap;
public:
    explicit CompareByPriority(phmap::flat_hash_map<std::string, int>& ComparatorMap) :ComparatorMap(ComparatorMap) {};

    // Comparator function
    bool operator()(std::shared_ptr<VariantContext>& vc1, std::shared_ptr<VariantContext>& vc2){
        assert(ComparatorMap.find(vc1->getSource()) != ComparatorMap.end());
        assert(ComparatorMap.find(vc2->getSource()) != ComparatorMap.end());
        return ComparatorMap[vc1->getSource()] < ComparatorMap[vc2->getSource()];
    }
};


// TODO: validate this method. In GATK, another vector is constructed for sorted.
std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> GATKVariantContextUtils::sortVariantContextsByPriority(
        std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> unsortedVCs,
        std::vector<std::string> &priorityListOfVCs, GenotypeMergeType mergeOption) {
    if ( priorityListOfVCs.empty() || mergeOption == GenotypeMergeType::UNSORTED )
        return unsortedVCs;

    phmap::flat_hash_map<std::string, int> ComparatorMap;
    for(int i=0; i<priorityListOfVCs.size(); i++)
    {
        if(ComparatorMap.find(priorityListOfVCs[i]) == ComparatorMap.end())
        {
            ComparatorMap.insert(std::pair<std::string, int>(priorityListOfVCs[i], i));
        }
    }
    CompareByPriority comparator(ComparatorMap);
    std::sort(unsortedVCs->begin(), unsortedVCs->end(), comparator);
    return unsortedVCs;
}

std::shared_ptr<Allele>
GATKVariantContextUtils::determineReferenceAllele(std::vector<std::shared_ptr<VariantContext>> &VCs) {
    return determineReferenceAllele(VCs, nullptr);
}

std::shared_ptr<Allele>
GATKVariantContextUtils::determineReferenceAllele(std::vector<std::shared_ptr<VariantContext>> &VCs,  Locatable* loc) {
    std::shared_ptr<Allele> ref = nullptr;
    for(auto& vc : VCs)
    {
        if(contextMatchesLoc(vc, loc)){
            auto myRef = vc->getReference();
            if(ref == nullptr || ref->getLength() < myRef->getLength())
                ref = myRef;
            else if(ref->getLength() == myRef->getLength() && !(*ref == *myRef))
                throw "The provided variant file(s) have inconsistent references for the same position";
        }
    }
    return ref;
}

bool GATKVariantContextUtils::isNonSymbolicExtendableAllele(std::shared_ptr<Allele> &allele) {
    return !(allele->getIsReference() || allele->getIsSymbolic() || (*allele) == (*Allele::SPAN_DEL));
}

bool GATKVariantContextUtils::contextMatchesLoc(std::shared_ptr<VariantContext> &vc, Locatable *loc)
{
    return loc == nullptr || loc->getStart() == vc->getStart();
}

std::shared_ptr<AlleleMapper> GATKVariantContextUtils::resolveIncompatibleAlleles(std::shared_ptr<Allele> refAllele,
                                                                                  std::shared_ptr<VariantContext> vc,
                                                                                  std::vector<std::shared_ptr<Allele>> &alleles) {
    if((*refAllele) == (*vc->getReference()))
        return std::make_shared<AlleleMapper>(vc);
    else{
        auto map = createAlleleMapping(refAllele, vc, alleles);
        map->insert({vc->getReference(), refAllele});
        return std::make_shared<AlleleMapper>(map);
    }
}

std::shared_ptr<std::map<std::shared_ptr<Allele>, std::shared_ptr<Allele>>>
GATKVariantContextUtils::createAlleleMapping(std::shared_ptr<Allele> refAllele, std::shared_ptr<VariantContext> oneVc,
                                             phmap::flat_hash_set<std::shared_ptr<Allele>> &currentAlleles) {
    auto myRef = oneVc->getReference();
    assert(refAllele->getLength() > myRef->getLength());

    int myRefLength = myRef->getLength();
    int refAlleleLength = refAllele->getLength();
    auto refAlleleBases = refAllele->getBases();
    std::shared_ptr<uint8_t[]> extraBases(new uint8_t[refAlleleLength - myRefLength]);
    memcpy(extraBases.get(), refAlleleBases.get() + myRefLength, refAlleleLength - myRefLength);

    std::shared_ptr<std::map<std::shared_ptr<Allele>, std::shared_ptr<Allele>>> map(new std::map<std::shared_ptr<Allele>, std::shared_ptr<Allele>>);
    for(auto a : oneVc->getAlternateAlleles())
    {
        if(isNonSymbolicExtendableAllele(a))
        {
            auto extended = Allele::extend(a, extraBases, refAlleleLength - myRefLength);
            for(auto & b : currentAlleles)
            {
                if(*extended == *b)
                {
                    extended = b;
                }
            }
            map->insert({a, extended});
        } else if(*a == *(Allele::SPAN_DEL)) {
            map->insert({a,a});
        }
    }
    return map;
}

std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>>
GATKVariantContextUtils::createAlleleMapping(std::shared_ptr<Allele> refAllele, std::shared_ptr<VariantContext> oneVc, const std::vector<std::shared_ptr<Allele>> &currentAlleles)
{
    auto myRef = oneVc->getReference();
    assert(refAllele->getLength() > myRef->getLength());

    int myRefLength = myRef->getLength();
    int refAlleleLength = refAllele->getLength();
    auto refAlleleBases = refAllele->getBases();
    std::shared_ptr<uint8_t[]> extraBases(new uint8_t[refAlleleLength - myRefLength]);
    memcpy(extraBases.get(), refAlleleBases.get() + myRefLength, refAlleleLength - myRefLength);

    std::shared_ptr<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>> map(new phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>);
    for(auto a : oneVc->getAlternateAlleles())
    {
        if(isNonSymbolicExtendableAllele(a))
        {
            auto extended = Allele::extend(a, extraBases, refAlleleLength - myRefLength);
            for(auto & b : currentAlleles)
            {
                if(*extended == *b)
                {
                    extended = b;
                }
            }
            map->insert({a, extended});
        } else if(*a == *(Allele::SPAN_DEL)) {
            map->insert({a,a});
        }
    }
    return map;
}

void GATKVariantContextUtils::mergeGenotypes(GenoTypesContext &mergedGenotypes, std::shared_ptr<VariantContext> oneVC,
                                             std::shared_ptr<AlleleMapper> alleleMapping, bool uniquifySamples) {
    int size = oneVC->getGenotypes()->getSize();
    for(int i=0; i<size; i++)
    {
        std::shared_ptr<Genotype> g = oneVC->getGenotypes()->get(i);
        std::string name = mergedSampleName(oneVC->getSource(), g->getSampleName(), uniquifySamples);
        if(!mergedGenotypes.containsSample(name)){
            // only add if the name is new
            std::shared_ptr<Genotype> newG = g;

            if(uniquifySamples || alleleMapping->needsRemapping())
            {
                auto allels = alleleMapping->needsRemapping() ? *alleleMapping->remap(g->getAlleles()): g->getAlleles();
                newG = GenotypeBuilder(g, name, allels).make(); //---remember to free this pointer
            }
            mergedGenotypes.add(newG);
        }
    }
}

std::string GATKVariantContextUtils::mergedSampleName(std::string trackName, std::string sampleName, bool uniquify) {
    return uniquify ? sampleName + "." + trackName : sampleName;
}

std::shared_ptr<VariantContext>
GATKVariantContextUtils::trimAlleles(std::shared_ptr<VariantContext> inputVC, bool trimForward, bool trimReverse) {
    assert(inputVC != nullptr);
    if ( inputVC->getNAlleles() <= 1 || inputVC->isSNP() )
        return inputVC;

    // see whether we need to trim common reference base from all alleles
    int revTrim = trimReverse ? computeReverseClipping(inputVC->getAlleles(), inputVC->getReference()->getBases(), inputVC->getReference()->getBasesLength()) : 0;
    auto revTrimVC = trimAlleles(inputVC, -1, revTrim);
    int fwdTrim = trimForward ? computeForwardClipping(revTrimVC->getAlleles()) : -1;
    return trimAlleles(revTrimVC, fwdTrim, 0);
}

std::shared_ptr<VariantContext>
GATKVariantContextUtils::trimAlleles(std::shared_ptr<VariantContext> inputVC, int fwdTrimEnd, int revTrim) {
    if( fwdTrimEnd == -1 && revTrim == 0 ) // nothing to do, so just return inputVC unmodified
        return inputVC;

    auto alleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
    auto originalToTrimmedAlleleMap = std::make_shared<phmap::flat_hash_map<std::shared_ptr<Allele>, std::shared_ptr<Allele>, hash_Allele, equal_Allele>>();

    for(auto a : inputVC->getAlleles())
    {
        if(a->getIsSymbolic()){
            alleles->emplace_back(a);
            originalToTrimmedAlleleMap->insert({a,a});
        } else {
            // get bases for current allele and create a new one with trimmed bases
            int newLength = a->getLength() - revTrim - fwdTrimEnd -1;
            std::shared_ptr<uint8_t[]> newBases(new uint8_t[newLength]);
            std::copy(a->getBases().get() + fwdTrimEnd + 1, a->getBases().get() + a->getLength() - revTrim, newBases.get());
            auto trimmedAllele = std::make_shared<Allele>(newBases, newLength, a->getIsReference());
            alleles->emplace_back(trimmedAllele);
            originalToTrimmedAlleleMap->insert({a, trimmedAllele});
        }
    }

    // now we can recreate new genotypes with trimmed alleles
    AlleleMapper alleleMapper(originalToTrimmedAlleleMap);
    auto genotypes = updateGenotypesWithMappedAlleles(inputVC->getGenotypes(), alleleMapper);

    int start = inputVC->getStart() + (fwdTrimEnd + 1);
    VariantContextBuilder builder(inputVC);
    builder.setStart(start);
    builder.setStop(start + alleles->operator[](0)->getLength() - 1);
    builder.setAlleles(alleles);
    builder.setGenotypes(genotypes);
    return builder.make();
}

int
GATKVariantContextUtils::computeReverseClipping(std::vector<std::shared_ptr<Allele>>& unclippedAlleles,
                                                std::shared_ptr<uint8_t[]> ref, int refLength) {
    int clipping = 0;
    bool stillClipping = true;

    while ( stillClipping ) {
        for (auto a : unclippedAlleles ) {
            if ( a->getIsSymbolic() )
                continue;

            // we need to ensure that we don't reverse clip out all of the bases from an allele because we then will have the wrong
            // position set for the VariantContext (although it's okay to forward clip it all out, because the position will be fine).
            if ( a->getLength() - clipping == 0 )
                return clipping - 1;

            if ( a->getLength() - clipping <= 0 || a->getLength() == 0 ) {
                stillClipping = false;
            }
            else if ( refLength == clipping ) {
                return -1;
            }
            else if ( a->getBases()[a->getLength()-clipping-1] != ref[refLength-clipping-1] ) {
                stillClipping = false;
            }
        }
        if ( stillClipping )
            clipping++;
    }

    return clipping;
}

int GATKVariantContextUtils::computeForwardClipping(std::vector<std::shared_ptr<Allele>> &unclippedAlleles) {
    // cannot clip unless there's at least 1 alt allele
    if ( unclippedAlleles.size() <= 1 )
        return -1;

    // we cannot forward clip any set of alleles containing a symbolic allele
    int minAlleleLength = INT_MAX;
    for (auto a : unclippedAlleles ) {
        if ( a->getIsSymbolic() )
            return -1;
        minAlleleLength = std::min(minAlleleLength, a->getLength());
    }

    auto firstAlleleBases = unclippedAlleles[0]->getBases();
    int indexOflastSharedBase = -1;

    // the -1 to the stop is that we can never clip off the right most base
    for ( int i = 0; i < minAlleleLength - 1; i++) {
        uint8_t base = firstAlleleBases[i];

        for (  auto allele : unclippedAlleles ) {
            if ( allele->getBases()[i] != base )
                return indexOflastSharedBase;
        }

        indexOflastSharedBase = i;
    }

    return indexOflastSharedBase;
}

std::shared_ptr<GenoTypesContext> GATKVariantContextUtils::updateGenotypesWithMappedAlleles(std::shared_ptr<GenoTypesContext> originalGenotypes,
                                                                            AlleleMapper &alleleMapper) {
    auto updatedGenotypes = std::make_shared<GenoTypesContext>(originalGenotypes->getSize());

    for(int i=0; i<originalGenotypes->getSize(); i++)
    {
        auto genotype = originalGenotypes->get(i);
        auto updatedAlleles = alleleMapper.remap(genotype->getAlleles());
        updatedGenotypes->add(GenotypeBuilder(genotype).setAlleles(updatedAlleles).make());
    }

    return updatedGenotypes;
}

std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>>
GATKVariantContextUtils::getNumTandemRepeatUnits(const std::shared_ptr<VariantContext> &vc,
                                                 const std::shared_ptr<uint8_t[]> refBasesStartingAtVCWithPad, int len) {
    if(! vc->isIndel()) {
        return {};
    }
    bool VERBOSE = false;
    int newLen;
    auto refBasesStartingAtVCWithoutPad = Mutect2Utils::copyOfRange(refBasesStartingAtVCWithPad, len, 1, len, newLen);
    auto refAllele = vc->getReference();
    int refLen;
    auto refAlleleBases = Mutect2Utils::copyOfRange(refAllele->getBases(), refAllele->getBasesLength(), 1, refAllele->getBasesLength(), refLen);
    std::shared_ptr<uint8_t[]> repeatUnit;
    std::vector<int> lengths;
    for(auto & allele : vc->getAlternateAlleles()) {
        int alleleLen;
        auto alleleBases = Mutect2Utils::copyOfRange(allele->getBases(), allele->getBasesLength(), 1, allele->getBasesLength(), alleleLen);
        auto result = getNumTandemRepeatUnits(refAlleleBases, refLen, alleleBases, alleleLen, refBasesStartingAtVCWithoutPad, newLen);
        std::vector<int> repetitionCount = result.first;
        if (repetitionCount[0] == 0 || repetitionCount[1] == 0)
            return {};

        if (lengths.empty()) {
            lengths.emplace_back(repetitionCount[0]); // add ref allele length only once
        }
        lengths.emplace_back(repetitionCount[1]);  // add this alt allele's length

        repeatUnit = result.second;
    }
    return {lengths, repeatUnit};
}

std::pair<std::vector<int>, std::shared_ptr<uint8_t[]>>
GATKVariantContextUtils::getNumTandemRepeatUnits(const std::shared_ptr<uint8_t[]> &refBases, int refLen,
                                                 const std::shared_ptr<uint8_t[]> &altBases, int altLen,
                                                 const std::shared_ptr<uint8_t[]> &remainingRefContext, int remainLen) {
    uint8_t * longB;
    int len;
    if(altLen > refLen) {
        longB = altBases.get();
        len = altLen;
    } else {
        longB = refBases.get();
        len = refLen;
    }

    int repeatUnitLength = findRepeatedSubstring(longB, refLen);
    std::shared_ptr<uint8_t[]> repeatUnit = std::shared_ptr<uint8_t[]>(new uint8_t[repeatUnitLength]);
    memcpy(repeatUnit.get(), longB, repeatUnitLength);
    std::vector<int> repetitionCount(2);
    int repetitionsInRef = findNumberOfRepetitions(repeatUnit.get(), repeatUnitLength, refBases.get(), refLen, true);
    uint8_t * tmp1 = new uint8_t [refLen+remainLen];
    memcpy(tmp1, refBases.get(), refLen);
    memcpy(tmp1+refLen, remainingRefContext.get(), remainLen);
    uint8_t * tmp2 = new uint8_t [altLen+remainLen];
    memcpy(tmp2, altBases.get(), altLen);
    memcpy(tmp2+altLen, remainingRefContext.get(), remainLen);
    repetitionCount[0] = findNumberOfRepetitions(repeatUnit.get(), repeatUnitLength, tmp1, refLen+remainLen,true)-repetitionsInRef;
    repetitionCount[1] = findNumberOfRepetitions(repeatUnit.get(), repeatUnitLength, tmp2, altLen+remainLen, true)-repetitionsInRef;
    delete[] tmp1;
    delete[] tmp2;
    return {repetitionCount, repeatUnit};
}

int GATKVariantContextUtils::findRepeatedSubstring(uint8_t *bases, int basesLen) {
    int repLength;
    for (repLength=1; repLength <= basesLen; repLength++) {
        uint8_t * candidateRepeatUnit = new uint8_t[repLength];
        memcpy(candidateRepeatUnit, bases, repLength);
        bool allBasesMatch = true;
        for (int start = repLength; start < basesLen; start += repLength ) {
            // check that remaining of string is exactly equal to repeat unit
            uint8_t * basePiece = new uint8_t[repLength];
            memcpy(basePiece, bases+start, repLength);
            if (!memcpy(candidateRepeatUnit, basePiece, repLength)) {
                allBasesMatch = false;
                break;
            }
            delete[] basePiece;
        }
        delete[] candidateRepeatUnit;
        if (allBasesMatch)
            return repLength;
    }

    return repLength;
}
