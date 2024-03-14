//
// Created by lhh on 5/20/22.
//

#ifndef MUTECT2CPP_MASTER_CALLEDHAPLOTYPES_H
#define MUTECT2CPP_MASTER_CALLEDHAPLOTYPES_H

#include "parallel_hashmap/phmap.h"
#include "Haplotype.h"
#include "VariantContext.h"

/**
 * Carries the result of a call to #assignGenotypeLikelihoods
 */
class CalledHaplotypes {
private:
    std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> calls;
    std::shared_ptr<phmap::flat_hash_set<std::shared_ptr<Haplotype>>> calledHaplotypes;

public:
    CalledHaplotypes(std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> calls, std::shared_ptr<phmap::flat_hash_set<std::shared_ptr<Haplotype>>> calledHaplotypes);

    /**
     * Get the list of calls made at this location
     * @return a non-null (but potentially empty) list of calls
     */
    std::shared_ptr<std::vector<std::shared_ptr<VariantContext>>> getCalls();

    /**
     * Get the set of haplotypes that we actually called (i.e., underlying one of the VCs in getCalls().
     * @return a non-null set of haplotypes
     */
    std::shared_ptr<phmap::flat_hash_set<std::shared_ptr<Haplotype>>> getCalledHaplotypes();
};


#endif //MUTECT2CPP_MASTER_CALLEDHAPLOTYPES_H
