//
// Created by lhh on 6/2/22.
//

#ifndef MUTECT2CPP_MASTER_PERALLELEANNOTATION_H
#define MUTECT2CPP_MASTER_PERALLELEANNOTATION_H

#include "InfoFieldAnnotation.h"

/**
 * Apply an annotation based on aggregation data from all reads supporting each allele.
 */
class PerAlleleAnnotation : public InfoFieldAnnotation{
private:
    bool  includeAllele(shared_ptr<Allele> allele);

protected:
    virtual std::optional<int> getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc) = 0;

    virtual int aggregate(std::vector<int>& values) = 0;

    virtual std::string getVcfKey() = 0;

    virtual std::string getDescription() = 0;

    // this is false by default but implementations may wish to override
    virtual bool includeRefAllele() { return false; }


public:
    /**
    * Calculate annotations for each allele based on given VariantContext and likelihoods for a given genotype's sample
    * and add the annotations to the GenotypeBuilder.  By default annotations are only calculated for alt alleles but
    * implementations may override the {@code includeRefAllele()} method.  See parent class docs in {@link GenotypeAnnotation}.
    */
    std::shared_ptr<std::map<std::string, AttributeValue>> annotate(shared_ptr<ReferenceContext> ref, shared_ptr<VariantContext> vc, AlleleLikelihoods<SAMRecord, Allele>* likelihoods);

    bool isUsableRead(std::shared_ptr<SAMRecord> read);

    std::vector<std::string> getKeyNames();


};


#endif //MUTECT2CPP_MASTER_PERALLELEANNOTATION_H
