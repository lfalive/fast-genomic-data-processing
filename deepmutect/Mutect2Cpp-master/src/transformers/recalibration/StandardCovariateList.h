#ifndef STANDARDCOVARIATE_LIST_H
#define STANDARDCOVARIATE_LIST_H

#include <vector>
#include "ContextCovariate.h"
#include "CycleCovariate.h"
#include "htslib/sam.h"
#include "SamRead.h"
#include "RecalibrationArgumentCollection.h"
using namespace std;

class StandardCovariateList{
private:
    vector<Covariate *> additionalCovariates;
    sam_hdr_t * readsHeader;
    Covariate * pContextCovariate;
    Covariate * pCycleCovariate;


public:
    const static int numberOfSpecialCovariate = 2;

    StandardCovariateList(sam_hdr_t * hdr, BaseArgument baseArgument);

    StandardCovariateList(RecalibrationArgumentCollection & RAC);

    ~StandardCovariateList();

    void recordAllValuesInStorage(bam1_t * read, Key * keys);

    int numberOfSpecialCovariates();

    int getMaximumKeyValue(int index);

    string getCovaariateName(int index);

    vector<Covariate *> & getAdditionalCovariates();

    /**
     * Retrieves a covariate by the parsed name or null
     * if no covariate with that name exists in the list.
     * Covariate name can only be "Context" or "Cycle"
     */
    Covariate * getCovariateByParsedName(string & covName);
};

#endif
