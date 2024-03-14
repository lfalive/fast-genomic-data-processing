#include <iostream>
#include <cassert>
#include "StandardCovariateList.h"

StandardCovariateList::StandardCovariateList(sam_hdr_t  *hdr, BaseArgument baseArgument) {
    this->readsHeader = hdr;
    this->pContextCovariate = new ContextCovariate(baseArgument);
    this->pCycleCovariate = new CycleCovariate(baseArgument);
    this->additionalCovariates.push_back(pContextCovariate);
    this->additionalCovariates.push_back(pCycleCovariate);
}

StandardCovariateList::StandardCovariateList(RecalibrationArgumentCollection &RAC)
{
    this->readsHeader = nullptr;    // unused variable
    this->pContextCovariate = new ContextCovariate(RAC);
    this->pCycleCovariate = new CycleCovariate(RAC);
    additionalCovariates.push_back(pContextCovariate);
    additionalCovariates.push_back(pCycleCovariate);
}

StandardCovariateList::~StandardCovariateList() {
    delete pContextCovariate;
    delete pCycleCovariate;
}

void StandardCovariateList::recordAllValuesInStorage( bam1_t *read, Key * keys){
    for (unsigned i=0; i<this->additionalCovariates.size(); i++)
    {
        this->additionalCovariates[i]->recordValues(read, readsHeader, keys);
    }

}

int StandardCovariateList::numberOfSpecialCovariates()
{
    return 2;
}

int StandardCovariateList::getMaximumKeyValue(int index)
{
    return additionalCovariates[index]->maximumKeyValue();
}

vector<Covariate *> & StandardCovariateList::getAdditionalCovariates()
{
    return additionalCovariates;
}

string StandardCovariateList::getCovaariateName(int index)
{
    assert(index >= 0 && index < additionalCovariates.size());
    return additionalCovariates[index]->getName();
}

Covariate *StandardCovariateList::getCovariateByParsedName(string &covName)
{
    if(covName == "Context")
        return pContextCovariate;
    else if(covName == "Cycle")
        return pCycleCovariate;
    else
        return nullptr;
}