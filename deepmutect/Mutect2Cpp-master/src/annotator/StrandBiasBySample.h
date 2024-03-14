//
// Created by lhh on 6/7/22.
//

#ifndef MUTECT2CPP_MASTER_STRANDBIASBYSAMPLE_H
#define MUTECT2CPP_MASTER_STRANDBIASBYSAMPLE_H


#include "GenotypeAnnotation.h"

class StrandBiasBySample : public GenotypeAnnotation {
    //For now this is only for 2x2 contingency tables
private:
    const static int ARRAY_DIM = 2;

public:
    void annotate(ReferenceContext &ref, shared_ptr<VariantContext> vc, std::shared_ptr<Genotype> g, GenotypeBuilder &gb, AlleleLikelihoods<SAMRecord, Allele> *likelihoods) override;

    /**
    * Helper function to turn the FisherStrand 2x2 table into the SB annotation array
    * @param table the 2x2 table used by the FisherStrand annotation
    * @return the array used by the per-sample Strand Bias annotation
    */
    static vector<int> getContingencyArray(vector<vector<int>>& table);

    std::vector<std::string> getKeyNames();
};


#endif //MUTECT2CPP_MASTER_STRANDBIASBYSAMPLE_H
