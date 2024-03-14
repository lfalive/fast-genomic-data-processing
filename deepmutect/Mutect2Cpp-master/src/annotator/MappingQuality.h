//
// Created by lhh on 6/15/22.
//

#ifndef MUTECT2CPP_MASTER_MAPPINGQUALITY_H
#define MUTECT2CPP_MASTER_MAPPINGQUALITY_H

#include "PerAlleleAnnotation.h"

class MappingQuality : public PerAlleleAnnotation {
private:
    const static int VALUE_FOR_NO_READS = 60;

protected:
    int aggregate(std::vector<int>& values);

    std::optional<int> getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc);

    std::string getVcfKey();

    std::string getDescription();

    // this is false by default but implementations may wish to override
    bool includeRefAllele() { return true; }

};


#endif //MUTECT2CPP_MASTER_MAPPINGQUALITY_H
