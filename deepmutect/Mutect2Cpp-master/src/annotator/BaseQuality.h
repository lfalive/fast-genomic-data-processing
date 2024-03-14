//
// Created by lhh on 6/2/22.
//

#ifndef MUTECT2CPP_MASTER_BASEQUALITY_H
#define MUTECT2CPP_MASTER_BASEQUALITY_H

#include "PerAlleleAnnotation.h"

/**
 * Median base quality of bases supporting each allele.
 *
 * <p>The output is an array containing, for each allele, the median base quality at the variant (for SNVs) and one base before the variant (for indels) over all reads that best match that allele.</p>
 *
 * <p>For example, for variant context with ref allele A and alt allele C we use base qualities at the A and C.  For variant context with ref allele AG and alt allele A (deletion),
 * we use base qualities at the A.  For variant context with ref allele A and alt allele AG (insertion) we use base qualities at the A.</p>
 */
class BaseQuality : public PerAlleleAnnotation{
protected:
    int aggregate(std::vector<int>& values);

    std::string getVcfKey();

    bool includeRefAllele() override;

    std::optional<int> getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc);

    std::string getDescription();
public:
    // -1 represent null
    std::optional<int> getBaseQuality(const std::shared_ptr<SAMRecord>& read, const shared_ptr<VariantContext>& vc);
};


#endif //MUTECT2CPP_MASTER_BASEQUALITY_H
