//
// Created by 梦想家xixi on 2021/12/8.
//

#include <limits>
#include <memory>
#include <builder/GenotypeBuilder.h>
#include "AssemblyBasedCallerUtils.h"
#include "haplotypecaller/ReferenceConfidenceModel.h"
#include "clipping/ReadClipper.h"
#include "read/ReadUtils.h"
#include "QualityUtils.h"
#include "utils/fragments/FragmentCollection.h"
#include "utils/fragments/FragmentUtils.h"
#include "AlignmentUtils.h"
#include "LocationAndAlleles.h"
#include "utils/variant/GATKVariantContextUtils.h"
#include "variantcontext/builder/VariantContextBuilder.h"

std::string AssemblyBasedCallerUtils::phase01 = "0|1";
std::string AssemblyBasedCallerUtils::phase10 = "1|0";

std::shared_ptr<Haplotype>
AssemblyBasedCallerUtils::createReferenceHaplotype(const std::shared_ptr<AssemblyRegion> &region,
                                                   const std::shared_ptr<SimpleInterval> &referencePadding,
                                                   ReferenceCache &cache) {
	int length = 0;
	std::shared_ptr<uint8_t[]> tmp = region->getAssemblyRegionReference(&cache, 0, length);
	return ReferenceConfidenceModel::createReferenceHaplotype(region, tmp, length, referencePadding);
}

std::shared_ptr<AssemblyResultSet>
AssemblyBasedCallerUtils::assembleReads(const std::shared_ptr<AssemblyRegion> &region,
                                        M2ArgumentCollection &argumentCollection,
                                        SAMFileHeader *header, ReferenceCache &cache,
                                        ReadThreadingAssembler &assemblyEngine) {
	finalizeRegion(region, false, false, 9, header, false);
	int refLength = 0;
	std::shared_ptr<uint8_t[]> fullReferenceWithPadding = region->getAssemblyRegionReference(&cache,
	                                                                                         REFERENCE_PADDING_FOR_ASSEMBLY,
	                                                                                         refLength);
	const std::shared_ptr<SimpleInterval> paddedReferenceLoc = getPaddedReferenceLoc(region,
	                                                                                 REFERENCE_PADDING_FOR_ASSEMBLY,
	                                                                                 header);
	std::shared_ptr<Haplotype> refHaplotype = createReferenceHaplotype(region, paddedReferenceLoc, cache);
	std::shared_ptr<AssemblyResultSet> assemblyResultSet = assemblyEngine.runLocalAssembly(region, refHaplotype,
	                                                                                       fullReferenceWithPadding,
	                                                                                       refLength,
	                                                                                       paddedReferenceLoc,
	                                                                                       nullptr);
	return assemblyResultSet;
}

std::shared_ptr<SimpleInterval>
AssemblyBasedCallerUtils::getPaddedReferenceLoc(const std::shared_ptr<AssemblyRegion> &region, int referencePadding,
                                                SAMFileHeader *header) {
	int padLeft = std::max(region->getExtendedSpan()->getStart() - referencePadding, 0);
	int padRight = std::min(region->getExtendedSpan()->getEnd() + referencePadding,
	                        header->getSequenceDictionary().getSequence(
			                        region->getExtendedSpan()->getContig()).getSequenceLength() - 1);
	return std::make_shared<SimpleInterval>(region->getExtendedSpan()->getContigInt(), padLeft, padRight);
}

void
AssemblyBasedCallerUtils::finalizeRegion(const std::shared_ptr<AssemblyRegion> &region, bool errorCorrectReads,
                                         bool dontUseSoftClippedBases,
                                         uint8_t minTailQuality, SAMFileHeader *header,
                                         bool correctOverlappingBaseQualities) {
	if (region->isFinalized())
		return;
	std::vector<std::shared_ptr<SAMRecord>> readsToUse;
	for (const std::shared_ptr<SAMRecord> &myRead: region->getReads()) {
		uint8_t minTailQualityToUse = errorCorrectReads ? 6 : minTailQuality;
		std::shared_ptr<SAMRecord> clippedRead = ReadClipper::hardClipLowQualEnds(myRead, minTailQualityToUse);
		clippedRead = dontUseSoftClippedBases || !ReadUtils::hasWellDefinedFragmentSize(clippedRead) ?
		              ReadClipper::hardClipSoftClippedBases(clippedRead) :
		              ReadClipper::revertSoftClippedBases(clippedRead);

		clippedRead = clippedRead->isUnmapped() ? clippedRead : ReadClipper::hardClipAdaptorSequence(clippedRead);
		if (!clippedRead->isEmpty() && clippedRead->getCigar()->getReadLength() > 0) {
			clippedRead = ReadClipper::hardClipToRegion(clippedRead, region->getExtendedSpan()->getStart(),
			                                            region->getExtendedSpan()->getEnd());
			if (region->readOverlapsRegion(clippedRead) && clippedRead->getLength() > 0) {
				readsToUse.emplace_back(
						(clippedRead == myRead) ? std::make_shared<SAMRecord>(*clippedRead)
						                        : clippedRead);
			}
		}
	}
	region->clearReads();
	region->addAll(readsToUse);
#ifdef SORT_MODE
	region->sortReadsByCoordinate();
#endif
	region->setFinalized(true);
}

std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>
AssemblyBasedCallerUtils::splitReadsBySample(const std::vector<std::string>& sampleList, const std::string& normalSample, const std::vector<std::shared_ptr<SAMRecord>> &reads) {
	std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> res = std::make_shared<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>();
	string tumorSample;
	for(auto& sample : sampleList)
	{
	    if(sample != normalSample)
	        tumorSample = sample;
	}

	res->insert({normalSample, std::vector<std::shared_ptr<SAMRecord>>()});
	res->insert({tumorSample, std::vector<std::shared_ptr<SAMRecord>>()});
	std::vector<std::shared_ptr<SAMRecord>> &normalReads = res->at(normalSample);
	std::vector<std::shared_ptr<SAMRecord>> &tumorReads = res->at(tumorSample);
	for (const std::shared_ptr<SAMRecord> &read: reads) {
		if (read->getGroup() == 0) {
			normalReads.emplace_back(read);
		} else {
			tumorReads.emplace_back(read);
		}
	}
	return res;
}

PairHMMLikelihoodCalculationEngine *
AssemblyBasedCallerUtils::createLikelihoodCalculationEngine(LikelihoodEngineArgumentCollection &likelihoodArgs) {
	double log10GlobalReadMismappingRate =
			likelihoodArgs.phredScaledGlobalReadMismappingRate < 0 ? (-1) * std::numeric_limits<double>::infinity()
			                                                       : QualityUtils::qualToErrorProbLog10(
					likelihoodArgs.phredScaledGlobalReadMismappingRate);
	return new PairHMMLikelihoodCalculationEngine((char) likelihoodArgs.gcpHMM, likelihoodArgs.pairHMMNativeArgs,
	                                              log10GlobalReadMismappingRate, likelihoodArgs.pcrErrorModel,
	                                              likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD);
}

void AssemblyBasedCallerUtils::cleanOverlappingReadPairs(vector<shared_ptr<SAMRecord>> &reads, const std::vector<std::string>& sampleList, const string &sample,
                                                         bool setConflictingToZero, int halfOfPcrSnvQual,
                                                         int halfOfPcrIndelQual) {
	auto MappedReads = splitReadsBySample(sampleList, sample, reads);
	for (auto &iter: *MappedReads) {
		FragmentCollection<SAMRecord> *fragmentCollection = FragmentCollection<SAMRecord>::create(iter.second);
		for (auto &overlappingPair: fragmentCollection->getOverlappingPairs()) {
			FragmentUtils::adjustQualsOfOverlappingPairedFragments(overlappingPair, setConflictingToZero,
			                                                       halfOfPcrSnvQual, halfOfPcrIndelQual);
		}
		delete fragmentCollection;
	}
}

std::shared_ptr<AssemblyRegion> AssemblyBasedCallerUtils::assemblyRegionWithWellMappedReads(
		const std::shared_ptr<AssemblyRegion> &originalAssemblyRegion, int minMappingQuality, SAMFileHeader *header) {
	auto result = make_shared<AssemblyRegion>(*originalAssemblyRegion->getSpan(),
	                                          originalAssemblyRegion->getSupportingStates(),
	                                          originalAssemblyRegion->getIsActive(),
	                                          originalAssemblyRegion->getExtension(), header);
	for (auto &read: originalAssemblyRegion->getReads()) {
		if (read->getMappingQuality() >= minMappingQuality)
			result->add(read);
	}
#ifdef SORT_MODE
	result->sortReadsByCoordinate();
#endif
	return result;
}

shared_ptr<phmap::flat_hash_map<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>> AssemblyBasedCallerUtils::realignReadsToTheirBestHaplotype(
        AlleleLikelihoods<SAMRecord, Haplotype> &originalReadLikelihoods, shared_ptr<Haplotype> &refHaplotype, shared_ptr<SimpleInterval>& paddedReferenceLoc,
        SmithWatermanAligner *aligner) {
    auto bestAlleles = originalReadLikelihoods.bestAllelesBreakingTies(&AssemblyBasedCallerUtils::HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY);
    auto result = std::make_shared<phmap::flat_hash_map<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>>(bestAlleles->size());

    for(auto & bestAllele : *bestAlleles)
    {
        auto & originalRead = bestAllele->evidence;
        auto & bestHaplotype = bestAllele->allele;
        bool isInformative = bestAllele->isInformative();
        shared_ptr<SAMRecord> realignedRead = AlignmentUtils::createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc->getStart(), isInformative, aligner);
        result->insert(pair<shared_ptr<SAMRecord>, shared_ptr<SAMRecord>>(originalRead, realignedRead));
    }
    return result;
}

double AssemblyBasedCallerUtils::HAPLOTYPE_ALIGNMENT_TIEBREAKING_PRIORITY(shared_ptr<Haplotype> h) {
    assert(h != nullptr);
    auto cigar = h->getCigar();
    int referenceTerm = h->getIsReference() ? 1 : 0;
    int cigarTerm = cigar == nullptr ? 0 : 1-cigar->numCigarElements();

    return (double)(referenceTerm + cigarTerm);
}

shared_ptr<vector<shared_ptr<VariantContext>>>
AssemblyBasedCallerUtils::getVariantContextsFromActiveHaplotypes(int loc, vector<shared_ptr<Haplotype>> &haplotypes,
                                                                 bool includeSpanningEvents) {
    shared_ptr<vector<shared_ptr<VariantContext>>> results = make_shared<vector<shared_ptr<VariantContext>>>();
    phmap::flat_hash_set<shared_ptr<LocationAndAlleles>, hash_LocationAndAlleles, equal_LocationAndAlleles> uniqueLocationsAndAlleles;

    // transform a haplotype to a stream of VariantContext

    vector<shared_ptr<VariantContext>> temp;
    for(auto& h : haplotypes)
    {
        auto events = h->getEventMap()->getOverlappingEvents(loc);
        for(shared_ptr<VariantContext>& event : *events)
        {
            temp.emplace_back(event);
        }
    }

    for(auto& v : temp)
    {
        if(v->getStart() == loc)
        {
            shared_ptr<LocationAndAlleles> locationAndAlleles = make_shared<LocationAndAlleles>(v->getStart(), v->getAlleles());
            if(uniqueLocationsAndAlleles.find(locationAndAlleles) == uniqueLocationsAndAlleles.end())
            {
                uniqueLocationsAndAlleles.insert(locationAndAlleles);
                results->emplace_back(v);
            }
        }
    }
    return results;
}

shared_ptr<VariantContext>
AssemblyBasedCallerUtils::makeMergedVariantContext(shared_ptr<vector<shared_ptr<VariantContext>>> vcs) {
    if(vcs->empty())
        return nullptr;

    vector<string> haplotypeSources;
    for(auto& variant : *vcs)
    {
        haplotypeSources.emplace_back(variant->getSource());
    }
    return GATKVariantContextUtils::simpleMerge(vcs, haplotypeSources, FilteredRecordMergeType::KEEP_IF_ANY_UNFILTERED, GenotypeMergeType::PRIORITIZE, false);
}

shared_ptr<phmap::flat_hash_map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Haplotype>>>, hash_Allele, equal_Allele>>
AssemblyBasedCallerUtils::createAlleleMapper(shared_ptr<VariantContext> mergedVC, int loc,
                                             vector<shared_ptr<Haplotype>> &haplotypes) {
    auto result = make_shared<phmap::flat_hash_map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Haplotype>>>, hash_Allele, equal_Allele>>();

    auto ref = mergedVC->getReference();
    result->insert({ref, make_shared<vector<shared_ptr<Haplotype>>>()});

    //Note: we can't use the alleles implied by eventsAtThisLoc because they are not yet merged to a common reference
    //For example, a homopolymer AAAAA reference with a single and double deletion would yield (i) AA* A and (ii) AAA* A
    //in eventsAtThisLoc, when in mergedVC it would yield AAA* AA A
    auto AlternateAlleles = mergedVC->getAlternateAlleles();
    for(auto& allele : AlternateAlleles)
    {
        if(!allele->getIsSymbolic())
            result->insert({allele, make_shared<vector<shared_ptr<Haplotype>>>()});
    }

    for(auto& h : haplotypes)
    {
        auto spanningEvents = h->getEventMap()->getOverlappingEvents(loc);
        if(spanningEvents->empty()){    //no events --> this haplotype supports the reference at this locus
            result->at(ref)->emplace_back(h);
            continue;
        }

        for(auto& spanningEvent : *spanningEvents)
        {
            if (spanningEvent->getStart() == loc) {
                // the event starts at the current location

                if (spanningEvent->getReference()->getLength() == mergedVC->getReference()->getLength()) {
                    // reference allele lengths are equal; we can just use the spanning event's alt allele
                    // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                    if(result->find(spanningEvent->getAlternateAllele(0)) != result->end())
                        result->at(spanningEvent->getAlternateAllele(0))->emplace_back(h);
                } else if(spanningEvent->getReference()->getLength() < mergedVC->getReference()->getLength())
                {
                    // spanning event has shorter ref allele than merged VC; we need to pad out its alt allele
                    phmap::flat_hash_set<std::shared_ptr<Allele>> currentAlleles;
                    auto spanningEventAlleleMappingToMergedVc = GATKVariantContextUtils::createAlleleMapping(mergedVC->getReference(), spanningEvent, currentAlleles);
                    auto remappedSpanningEventAltAllele = spanningEventAlleleMappingToMergedVc->at(spanningEvent->getAlternateAllele(0));
                    // in the case of GGA mode the spanning event might not match an allele in the mergedVC
                    if(result->find(remappedSpanningEventAltAllele) != result->end())
                    {
                        result->at(remappedSpanningEventAltAllele)->emplace_back(h);
                    }
                } else {
                    // the process of creating the merged VC in AssemblyBasedCallerUtils::makeMergedVariantContext should have
                    // already padded out the reference allele, therefore this spanning VC must not be in events at this site
                    // because we're in GGA mode and it's not an allele we want
                    continue;
                }
            } else {
                // the event starts prior to the current location, so it's a spanning deletion
                if(result->find(Allele::SPAN_DEL) == result->end())
                    result->insert({Allele::SPAN_DEL, make_shared<vector<shared_ptr<Haplotype>>>()});
                result->at(Allele::SPAN_DEL)->emplace_back(h);
                break;
            }
        }
    }
    return result;

}

shared_ptr<vector<shared_ptr<VariantContext>>>
AssemblyBasedCallerUtils::phaseCalls(vector<shared_ptr<VariantContext>> &calls,
                                     phmap::flat_hash_set<shared_ptr<Haplotype>> &calledHaplotypes) {
    // construct a mapping from alternate allele to the set of haplotypes that contain that allele
    auto haplotypeMap = constructHaplotypeMapping(calls, calledHaplotypes);

    map<VariantContext*, pair<int, string>> phaseSetMapping;
    int uniqueCounterEndValue = constructPhaseSetMapping(calls, *haplotypeMap, calledHaplotypes.size() - 1, phaseSetMapping);

    return constructPhaseGroups(calls, phaseSetMapping, uniqueCounterEndValue);
}

// ---a helper method
bool contains(vector<shared_ptr<Allele>> alleles, shared_ptr<Allele> alt)
{
    for(auto allele : alleles)
    {
        if(*allele == *alt)
            return true;
    }
    return false;
}

shared_ptr<map<VariantContext*, shared_ptr<phmap::flat_hash_set<Haplotype*>>>>
AssemblyBasedCallerUtils::constructHaplotypeMapping(vector<shared_ptr<VariantContext>> &originalCalls,
                                                    phmap::flat_hash_set<shared_ptr<Haplotype>> &calledHaplotypes) {
    auto haplotypeMap = make_shared<map<VariantContext*, shared_ptr<phmap::flat_hash_set<Haplotype*>>>>();
    for(auto call : originalCalls)
    {
        // don't try to phase if there is not exactly 1 alternate allele
        if(!isBiallelic(call)){
            haplotypeMap->insert({call.get(), nullptr});
            continue;
        }

        // keep track of the haplotypes that contain this particular alternate allele
        auto alt = call->getAlternateAllele(0);

        shared_ptr<phmap::flat_hash_set<Haplotype*>> hapsWithAllele = make_shared<phmap::flat_hash_set<Haplotype*>>();
        for(auto h : calledHaplotypes)
        {
            bool match = false;
            auto VariantContexts = h->getEventMap()->getVariantContexts();
            for(auto & varaint : VariantContexts)
            {
                if(varaint->getStart() == call->getStart() && contains(varaint->getAlternateAlleles(), alt))
                    match = true;
                if(*alt == *Allele::SPAN_DEL && varaint->getStart() < call->getStart() && varaint->getEnd() >= call->getStart())
                    match = true;
                if(match)
                    break;
            }
            if(match)
                hapsWithAllele->insert(h.get());
        }
        haplotypeMap->insert({call.get(), hapsWithAllele});
    }
    return haplotypeMap;
}

bool AssemblyBasedCallerUtils::isBiallelic(shared_ptr<VariantContext> vc) {
    if(vc->isBiallelic())
        return true;

    if(vc->getNAlleles() == 3)
    {
        for (auto& allele: vc->getAlternateAlleles()) {
            if(*allele == *Allele::NON_REF_ALLELE)
                return true;
        }
        return false;
    }
    return false;
}

//---whether set1 contain all of the elements of set2
bool containAll(phmap::flat_hash_set<Haplotype*>& set1, phmap::flat_hash_set<Haplotype*>& set2)
{
    for(auto h : set2)
    {
        if(set1.find(h) == set1.end())
            return false;
    }
    return true;
}

//---just like Java set retainAll
void retainAll(phmap::flat_hash_set<Haplotype*>& set1, phmap::flat_hash_set<Haplotype*>& set2)
{
    for (auto h = set1.begin(), last = set1.end(); h != last; ) {
        if (set2.find(*h) == set2.end()) {
            h = set1.erase(h);
        } else {
            ++h;
        }
    }
}

int AssemblyBasedCallerUtils::constructPhaseSetMapping(vector<shared_ptr<VariantContext>> &originalCalls,
                                                       map<VariantContext*, shared_ptr<phmap::flat_hash_set<Haplotype*>>> &haplotypeMap,
                                                       int totalAvailableHaplotypes,
                                                       map<VariantContext*, pair<int, string>> &phaseSetMapping) {
    int numCalls = originalCalls.size();
    int uniqueCounter = 0;

    // use the haplotype mapping to connect variants that are always/never present on the same haplotypes
    for ( int i = 0; i < numCalls - 1; i++ ) {
        shared_ptr<VariantContext>& call = originalCalls[i];
        shared_ptr<phmap::flat_hash_set<Haplotype*>> haplotypesWithCall = haplotypeMap.at(call.get());  // this variable can be nullptr
        if(haplotypesWithCall == nullptr || haplotypesWithCall->empty())
            continue;

        bool callIsOnAllHaps = haplotypesWithCall->size() == totalAvailableHaplotypes;
        for ( int j = i+1; j < numCalls; j++ ) {
            shared_ptr<VariantContext> comp = originalCalls[j];
            if(haplotypeMap.find(comp.get()) == haplotypeMap.end())
                continue;
            auto haplotypesWithComp = haplotypeMap.at(comp.get());
            if ( haplotypesWithComp == nullptr || haplotypesWithComp->empty()) {
                continue;
            }

            // if the variants are together on all haplotypes, record that fact.
            // another possibility is that one of the variants is on all possible haplotypes (i.e. it is homozygous).
            bool compIsOnAllHaps = haplotypesWithComp->size() == totalAvailableHaplotypes;
            if ( (haplotypesWithCall->size() == haplotypesWithComp->size() && containAll(*haplotypesWithCall ,*haplotypesWithComp)) || callIsOnAllHaps || compIsOnAllHaps ) {
                // create a new group if these are the first entries
                if ( phaseSetMapping.find(call.get()) ==  phaseSetMapping.end() ) {
                    // note that if the comp is already in the map then that is very bad because it means that there is
                    // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                    // situation, we should abort if we encounter it.
                    if ( phaseSetMapping.find(comp.get()) != phaseSetMapping.end() ) {
                        phaseSetMapping.clear();
                        return 0;
                    }

                    // An important note: even for homozygous variants we are setting the phase as "0|1" here.
                    // We do this because we cannot possibly know for sure at this time that the genotype for this
                    // sample will actually be homozygous downstream: there are steps in the pipeline that are liable
                    // to change the genotypes.  Because we can't make those assumptions here, we have decided to output
                    // the phase as if the call is heterozygous and then "fix" it downstream as needed.
                    phaseSetMapping.insert({call.get(), {uniqueCounter, phase01}});
                    phaseSetMapping.insert({comp.get(), {uniqueCounter, phase01}});
                    uniqueCounter++;
                }
                // otherwise it's part of an existing group so use that group's unique ID
                else if ( phaseSetMapping.find(comp.get()) ==  phaseSetMapping.end() ) {
                    pair<int, string> callPhase = phaseSetMapping.at(call.get());
                    phaseSetMapping.insert({comp.get(), {callPhase.first, callPhase.second}});
                }
            }
            // if the variants are apart on *all* haplotypes, record that fact
            else if ( haplotypesWithCall->size() + haplotypesWithComp->size() == totalAvailableHaplotypes ) {
                phmap::flat_hash_set<Haplotype*> intersection;
                intersection = *haplotypesWithCall;
                retainAll(intersection, *haplotypesWithComp);
                //intersection.retainAll(haplotypesWithComp);
                if ( intersection.empty() ) {
                    // create a new group if these are the first entries
                    if ( phaseSetMapping.find(call.get()) ==  phaseSetMapping.end() ) {
                        // note that if the comp is already in the map then that is very bad because it means that there is
                        // another variant that is in phase with the comp but not with the call.  Since that's an un-phasable
                        // situation, we should abort if we encounter it.
                        if (phaseSetMapping.find(comp.get()) != phaseSetMapping.end()) {
                            phaseSetMapping.clear();
                            return 0;
                        }

                        phaseSetMapping.insert({call.get(), {uniqueCounter, phase01}});
                        phaseSetMapping.insert({comp.get(), {uniqueCounter, phase10}});
                        uniqueCounter++;
                    }
                    // otherwise it's part of an existing group so use that group's unique ID
                    else if ( phaseSetMapping.find(comp.get()) == phaseSetMapping.end() ){
                        auto callPhase = phaseSetMapping.at(call.get());
                        phaseSetMapping.insert({comp.get(), {callPhase.first, callPhase.second == phase01 ? phase10 : phase01}});
                    }
                }
            }
        }
    }
    return uniqueCounter;

}

shared_ptr<vector<shared_ptr<VariantContext>>> AssemblyBasedCallerUtils::constructPhaseGroups(vector<shared_ptr<VariantContext>> &originalCalls, map<VariantContext*, pair<int, string>>& phaseSetMapping, int indexTo)
{
    auto phasedCalls = make_shared<vector<shared_ptr<VariantContext>>>(originalCalls);

    // if we managed to find any phased groups, update the VariantContexts
    for ( int count = 0; count < indexTo; count++ ) {
        // get all of the (indexes of the) calls that belong in this group (keeping them in the original order)
        vector<int> indexes;
        for ( int index = 0; index < originalCalls.size(); index++ ) {
            auto call = originalCalls[index];
            if ( phaseSetMapping.find(call.get()) != phaseSetMapping.end() && phaseSetMapping.at(call.get()).first == count ) {
                indexes.push_back(index);
            }
        }

        if ( indexes.size() < 2 ) {
            throw "Somehow we have a group of phased variants that has fewer than 2 members";
        }

        // create a unique ID based on the leftmost one
        string uniqueID = createUniqueID(originalCalls[indexes[0]]);

        // create the phase set identifier, which is the position of the first variant in the set
        int phaseSetID = originalCalls[indexes[0]]->getStart();

        // update the VCs
        for ( int index : indexes ) {
            auto originalCall = originalCalls[index];
            auto phasedCall = phaseVC(originalCall, uniqueID, phaseSetMapping.at(originalCall.get()).second, phaseSetID);
            phasedCalls->operator[](index) = phasedCall;
        }
    }
    return phasedCalls;
}

// TODO: Allele::getDisplayString equals to Allele::getBaseString() ?
string AssemblyBasedCallerUtils::createUniqueID(shared_ptr<VariantContext> vc) {
    return to_string(vc->getStart()) + "_" + vc->getReference()->getBaseString() + "_" + vc->getAlternateAllele(0)->getBaseString();
}

shared_ptr<VariantContext>
AssemblyBasedCallerUtils::phaseVC(shared_ptr<VariantContext> vc, string &ID, string &phaseGT, int phaseSetID) {
    vector<std::shared_ptr<Genotype>> phasedGenotypes;
    auto genotypesContext = vc->getGenotypes();
    if(genotypesContext != nullptr && genotypesContext->getGenotypes() != nullptr)
    {
        for(auto g : *genotypesContext->getGenotypes())
        {
            std::vector<std::shared_ptr<Allele>> & alleles = g->getAlleles();
            if (phaseGT == phase10 && g->isHet())
            {
                // reverse the list
                int size = alleles.size();
                for(int i=0; i<size/2; i++)
                {
                    swap(alleles[i], alleles[size - i - 1]);
                }
            }

            auto genotype = GenotypeBuilder(g).setAlleles(alleles).phased(true).attribute(VCFConstants::HAPLOTYPE_CALLER_PHASING_ID_KEY, ID).attribute(VCFConstants::HAPLOTYPE_CALLER_PHASING_GT_KEY, phaseGT).attribute(VCFConstants::PHASE_SET_KEY, phaseSetID).make();
            phasedGenotypes.push_back(genotype);
        }
    }
    return VariantContextBuilder(vc).setGenotypes(phasedGenotypes)->make();
}