//
// Created by lhh on 6/14/22.
//

#ifndef MUTECT2CPP_MASTER_FRAGMENTLENGTH_H
#define MUTECT2CPP_MASTER_FRAGMENTLENGTH_H

#include "PerAlleleAnnotation.h"

class FragmentLength : public PerAlleleAnnotation{
protected:
    int aggregate(std::vector<int>& values);

    bool includeRefAllele();

    std::string getVcfKey();

    std::string getDescription();

    std::optional<int> getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc);
};


#endif //MUTECT2CPP_MASTER_FRAGMENTLENGTH_H
