//
// Created by lhh on 5/24/22.
//

#ifndef MUTECT2CPP_MASTER_VARIANTCONTEXTUTILS_H
#define MUTECT2CPP_MASTER_VARIANTCONTEXTUTILS_H

#include "VariantContext.h"

class VariantContextUtils {
public:
    static int getSize(std::shared_ptr<VariantContext>& vc);

    /**
    * Update the attributes of the attributes map given the VariantContext to reflect the
    * proper chromosome-based VCF tags
    *
    * @param vc          the VariantContext
    * @param attributes  the attributes map to populate; must not be null; may contain old values
    * @param removeStaleValues should we remove stale values from the mapping?
    * @return the attributes map provided as input, returned for programming convenience
    */
    static void calculateChromosomeCounts(std::shared_ptr<VariantContext> vc, std::map<std::string, std::string>& attributes, bool removeStaleValues);
};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXTUTILS_H
