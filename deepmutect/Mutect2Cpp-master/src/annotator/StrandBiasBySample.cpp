//
// Created by lhh on 6/7/22.
//

#include "StrandBiasBySample.h"
#include "StrandBiasTest.h"

void
StrandBiasBySample::annotate(ReferenceContext &ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder &gb,
                             AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);
    assert(g != nullptr);

    if(likelihoods == nullptr || !g->isCalled())
    {
        std::cerr << "Annotation will not be calculated, genotype is not called or alleleLikelihoodMap is null\n";
        return;
    }

    auto table = StrandBiasTest::getContingencyTable(likelihoods, vc, 0, {g->getSampleName()});

    gb.attribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY, getContingencyArray(table));
}

vector<int> StrandBiasBySample::getContingencyArray(vector<vector<int>>& table)
{
    assert(table.size() == ARRAY_DIM && table[0].size() == ARRAY_DIM);
    vector<int> list;
    list.reserve(ARRAY_DIM * ARRAY_DIM);
    list.push_back(table[0][0]);
    list.push_back(table[0][1]);
    list.push_back(table[1][0]);
    list.push_back(table[1][1]);
    return list;
}

std::vector<std::string> StrandBiasBySample::getKeyNames() {
    return {VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY};
}