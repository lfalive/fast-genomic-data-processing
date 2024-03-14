//
// Created by lhh on 6/7/22.
//

#include "DepthPerSampleHC.h"

template<> double AlleleLikelihoods<SAMRecord, Allele>::NATURAL_LOG_INFORMATIVE_THRESHOLD = MathUtils::log10ToLog(LOG_10_INFORMATIVE_THRESHOLD);

void DepthPerSampleHC::annotate(ReferenceContext &ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder &gb,
                                AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);
    assert(g != nullptr);

    if(likelihoods == nullptr || !g->isCalled())
    {
        std::cerr << "Annotation will not be calculated, genotype is not called or alleleLikelihoodMap is null\n";
        return;
    }

    // check that there are reads
    string sample = g->getSampleName();
    if(likelihoods->sampleEvidenceCount(likelihoods->indexOfSamples(sample)) == 0){
        gb.setDP(0);
        return;
    }

    // the depth for the HC is the sum of the informative alleles at this site.  It's not perfect (as we cannot
    // differentiate between reads that align over the event but aren't informative vs. those that aren't even
    // close) but it's a pretty good proxy and it matches with the AD field (i.e., sum(AD) = DP).
    auto alleleSubset = make_shared<std::map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Allele>>>>>();
    for(auto allele: vc->getAlleles())
    {
        alleleSubset->insert({allele, make_shared<vector<shared_ptr<Allele>>>(1, allele)});
    }
    auto subsettedLikelihoods = likelihoods->marginalize(alleleSubset);
    auto bestAlleles = subsettedLikelihoods->bestAllelesBreakingTies(sample);
	delete subsettedLikelihoods;
    int depth = 0;

    for(auto ba : *bestAlleles)
    {
        if(ba->isInformative())
            depth++;
    }

    gb.setDP(depth);
}

std::vector<std::string> DepthPerSampleHC::getKeyNames() {
    return {VCFConstants::DEPTH_KEY};
}