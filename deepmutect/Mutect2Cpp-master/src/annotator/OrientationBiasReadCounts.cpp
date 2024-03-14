//
// Created by lhh on 6/7/22.
//

#include "read/ReadUtils.h"
#include "OrientationBiasReadCounts.h"

bool OrientationBiasReadCounts::isUsableRead(shared_ptr<SAMRecord> read) {
    return read->getMappingQuality() != 0 && read->getMappingQuality() != QualityUtils::MAPPING_QUALITY_UNAVALIABLE;
}

int OrientationBiasReadCounts::getReadBaseQuality(shared_ptr<SAMRecord> read, int refLoc) {
    assert(read != nullptr);
    int readCoord = ReadUtils::getReadCoordinateForReferenceCoordinate(read->getSoftStart(), read->getCigar(), refLoc, ClippingTail::RIGHT_TAIL, true);
    return readCoord < 0 || readCoord >= read->getLength() ? 0 : read->getBaseQuality(readCoord);
}

void OrientationBiasReadCounts::annotate(ReferenceContext &ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g,
                                         GenotypeBuilder &gb, AlleleLikelihoods<SAMRecord, Allele> *likelihoods) {
    assert(vc != nullptr);

    if(g == nullptr || likelihoods == nullptr)
        return;

    map<Allele*, int> f1r2Counts, f2r1Counts;
    for(auto allele : likelihoods->getAlleles())
    {
        f1r2Counts.insert({allele.get(), 0});
    }
    f2r1Counts = f1r2Counts;

    auto bestAlleles = likelihoods->bestAllelesBreakingTies(g->getSampleName());
    for(auto ba : *bestAlleles)
    {
        if(ba->isInformative() && isUsableRead(ba->evidence) && getReadBaseQuality(ba->evidence, vc->getStart()) >= MINIMUM_BASE_QUALITY)
        {
            if(ReadUtils::isF2R1(ba->evidence))
                f2r1Counts[ba->allele.get()]++;
            else
                f1r2Counts[ba->allele.get()]++;
        }
    }

    vector<int> f1r2;
    vector<int> f2r1;
    for(auto a : vc->getAlleles())
    {
        f1r2.push_back(f1r2Counts[a.get()]);
        f2r1.push_back(f2r1Counts[a.get()]);
    }

    gb.attribute(VCFConstants::F1R2_KEY, f1r2);
    gb.attribute(VCFConstants::F2R1_KEY, f2r1);
}

std::vector<std::string> OrientationBiasReadCounts::getKeyNames() {
    return {VCFConstants::F1R2_KEY, VCFConstants::F2R1_KEY};
}