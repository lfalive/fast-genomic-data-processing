//
// Created by lhh on 6/6/22.
//

#include "DepthPerAlleleBySample.h"

void
DepthPerAlleleBySample::annotate(ReferenceContext &ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder &gb,
                                 AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    if(g == nullptr || !g->isCalled() || likelihoods == nullptr)
        return;

    gb.setAD(annotateWithLikelihoods(vc, g, vc->getAlleles(), likelihoods));
}

std::vector<int> DepthPerAlleleBySample::annotateWithLikelihoods(shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g,
                                                                 vector<shared_ptr<Allele>> &alleles,
                                                                 AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    phmap::flat_hash_map<Allele*, int> alleleCounts;
    for(const auto& allele : vc->getAlleles())
    {
        alleleCounts.insert({allele.get(), 0});
    }

    auto alleleSubset = make_shared<std::map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Allele>>>>>();
    for(const auto& allele: alleles)
    {
        alleleSubset->insert({allele, make_shared<vector<shared_ptr<Allele>>>(1, allele)});
    }
    auto subsettedLikelihoods = likelihoods->marginalize(alleleSubset);
    auto bestAllels = subsettedLikelihoods->bestAllelesBreakingTies(g->getSampleName());
	delete subsettedLikelihoods;
	for(const auto& ba : *bestAllels)
    {
        if(ba->isInformative())
        {
            alleleCounts[ba->allele.get()]++;
        }
    }

	std::vector<int> counts(alleleCounts.size());
    counts[0] = alleleCounts[vc->getReference().get()];
    for (int i = 0; i < vc->getNAlleles() -1; i++) {
        counts[i + 1] = alleleCounts[vc->getAlternateAllele(i).get()];
    }

    return counts;
}

std::vector<std::string> DepthPerAlleleBySample::getKeyNames() {
    return {VCFConstants::GENOTYPE_ALLELE_DEPTHS};
}