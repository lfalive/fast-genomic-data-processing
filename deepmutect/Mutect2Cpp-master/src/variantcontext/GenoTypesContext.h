//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
#define MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
#include "parallel_hashmap/phmap.h"
#include "Genotype.h"

class Genotype;

class GenoTypesContext {
protected:
    std::vector<std::string> * sampleNamesInOrder;
    phmap::flat_hash_map<std::string, int> * sampleNameToOffset;
    std::vector<std::shared_ptr<Genotype>> * notToBeDirectlyAccessedGenotypes;


    void ensureSampleNameMap();
    void invalidateSampleOrdering();

private:
    int maxPloidy;
    bool immutable = false;
    std::vector<std::shared_ptr<Genotype>> bases;

public:
    static std::shared_ptr<GenoTypesContext> NO_GENOTYPES;

    explicit GenoTypesContext(int n = 10);
    explicit GenoTypesContext(std::vector<std::shared_ptr<Genotype>> & genotypes);
    GenoTypesContext(std::vector<std::shared_ptr<Genotype>> * genotypes, phmap::flat_hash_map<std::string, int> * sampleNameToOffset, std::vector<std::string> * sampleNamesInOrder);
    GenoTypesContext(std::vector<std::shared_ptr<Genotype>> * genotypes, phmap::flat_hash_map<std::string, int> * sampleNameToOffset, std::vector<std::string> * sampleNamesInOrder, bool immutable);
    ~GenoTypesContext();

    std::vector<std::shared_ptr<Genotype>> * getGenotypes();
    void setImmutable();
    int getSize();
    std::shared_ptr<Genotype> get(int i);
    bool containsSample(std::string sample);
    void add(std::shared_ptr<Genotype> genotype);
    bool isEmpty();
};


#endif //MUTECT2CPP_MASTER_GENOTYPESCONTEXT_H
