//
// Created by lhh on 11/12/21.
//

#include "GenoTypesContext.h"
#include <utility>
#include "parallel_hashmap/phmap.h"
#include <cassert>


GenoTypesContext::GenoTypesContext(int n) : maxPloidy(-1), immutable(false),  notToBeDirectlyAccessedGenotypes(new std::vector<std::shared_ptr<Genotype>>()), sampleNameToOffset(
        nullptr), sampleNamesInOrder(nullptr){
    notToBeDirectlyAccessedGenotypes->reserve(n);
}

GenoTypesContext::GenoTypesContext(std::vector<std::shared_ptr<Genotype>> & genotypes) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(
        new std::vector<std::shared_ptr<Genotype>>(genotypes)), sampleNameToOffset(nullptr), sampleNamesInOrder(nullptr){}

std::shared_ptr<GenoTypesContext> GenoTypesContext::NO_GENOTYPES = std::make_shared<GenoTypesContext>(nullptr, nullptr, nullptr, true);

GenoTypesContext::GenoTypesContext(std::vector<std::shared_ptr<Genotype>> *genotypes, phmap::flat_hash_map<std::string, int> *sampleNameToOffset,
                                   std::vector<std::string> *sampleNamesInOrder) : maxPloidy(-1), immutable(false), notToBeDirectlyAccessedGenotypes(genotypes), sampleNamesInOrder(sampleNamesInOrder), sampleNameToOffset(sampleNameToOffset){}

GenoTypesContext::GenoTypesContext(std::vector<std::shared_ptr<Genotype>> *genotypes, phmap::flat_hash_map<std::string, int> *sampleNameToOffset,
                                                                      std::vector<std::string> *sampleNamesInOrder, bool immutable) : maxPloidy(-1), immutable(immutable), notToBeDirectlyAccessedGenotypes(genotypes), sampleNamesInOrder(sampleNamesInOrder), sampleNameToOffset(sampleNameToOffset){}

GenoTypesContext::~GenoTypesContext() noexcept {
/*    if(notToBeDirectlyAccessedGenotypes)
    {
        for(Genotype* genotype : *notToBeDirectlyAccessedGenotypes)
        {
            delete genotype;
        }
    }*/

    delete notToBeDirectlyAccessedGenotypes;
    delete sampleNameToOffset;
    delete sampleNamesInOrder;
}

void GenoTypesContext::setImmutable() {
    immutable = true;
}

std::vector<std::shared_ptr<Genotype>> *GenoTypesContext::getGenotypes() {
    return notToBeDirectlyAccessedGenotypes;
}

int GenoTypesContext::getSize() {
    if(!notToBeDirectlyAccessedGenotypes)
        return 0;
    else
        return (int)notToBeDirectlyAccessedGenotypes->size();
}

std::shared_ptr<Genotype> GenoTypesContext::get(int i) {
    assert(notToBeDirectlyAccessedGenotypes != nullptr);
    assert(i >= 0 && i < notToBeDirectlyAccessedGenotypes->size());
    return (*notToBeDirectlyAccessedGenotypes)[i];
}

void GenoTypesContext::ensureSampleNameMap() {
    if(sampleNameToOffset == nullptr)
    {
        sampleNameToOffset = new phmap::flat_hash_map<std::string, int>(getSize());

        for(int i=0; i<getSize(); i++)
        {
            sampleNameToOffset->insert({notToBeDirectlyAccessedGenotypes->operator[](i)->getSampleName(), i});
        }
    }
}

bool GenoTypesContext::containsSample(std::string sample) {
    ensureSampleNameMap();
    return sampleNameToOffset->find(sample) != sampleNameToOffset->end();
}

void GenoTypesContext::add(std::shared_ptr<Genotype> genotype) {
    assert(!this->immutable);
    invalidateSampleOrdering();
    if(sampleNameToOffset != nullptr)
    {
        sampleNameToOffset->insert({genotype->getSampleName(), getSize()});
    }
    if(!notToBeDirectlyAccessedGenotypes)
        notToBeDirectlyAccessedGenotypes = new std::vector<std::shared_ptr<Genotype>>();
    notToBeDirectlyAccessedGenotypes->emplace_back(genotype);
}

void GenoTypesContext::invalidateSampleOrdering() {
    sampleNamesInOrder = nullptr;
}

bool GenoTypesContext::isEmpty() {
    return notToBeDirectlyAccessedGenotypes == nullptr || notToBeDirectlyAccessedGenotypes->empty();
}