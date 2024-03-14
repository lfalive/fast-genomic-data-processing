//
// Created by lhh on 6/15/22.
//

#ifndef MUTECT2CPP_MASTER_READPOSITION_H
#define MUTECT2CPP_MASTER_READPOSITION_H

#include "PerAlleleAnnotation.h"

class ReadPosition : public PerAlleleAnnotation{
private:
    const static int VALUE_FOR_NO_READS = 50;

protected:
    int aggregate(std::vector<int>& values);

    std::optional<int> getValueForRead(std::shared_ptr<SAMRecord> read, shared_ptr<VariantContext> vc);

    std::string getVcfKey();

    std::string getDescription();

};


#endif //MUTECT2CPP_MASTER_READPOSITION_H
