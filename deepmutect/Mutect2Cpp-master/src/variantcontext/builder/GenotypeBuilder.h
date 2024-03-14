//
// Created by 梦想家xixi on 2021/11/30.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPEBUILDER_H
#define MUTECT2CPP_MASTER_GENOTYPEBUILDER_H

#include <utility>
#include <vector>
#include "AttributeValue.h"
#include "FastGenotype.h"

class GenotypeBuilder {
private:
    static std::vector<std::shared_ptr<Allele>> HAPLOID_NO_CALL;
    static std::vector<std::shared_ptr<Allele>> DIPLOID_NO_CALL;
    std::string sampleName;
    std::vector<std::shared_ptr<Allele>> alleles;
    bool isPhased = false;
    int GQ = -1;
    int DP = -1;
	std::vector<int> AD;
	std::vector<int> PL;
    std::map<std::string, AttributeValue> extendedAttributes;
    std::string filters;
    int initialAttributeMapSize = 5;
    static std::map<std::string, AttributeValue> NO_ATTRIBUTES;

public:
    static std::shared_ptr<Genotype> create(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles);
    static std::shared_ptr<Genotype> create(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles, const std::map<std::string, AttributeValue>& attributes);
    static std::shared_ptr<Genotype> create(const std::string& sampleName, const std::vector<std::shared_ptr<Allele>>& alleles, std::vector<double> gls);


    GenotypeBuilder(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles) : sampleName(std::move(sampleName)), alleles(std::move(alleles)){}
    GenotypeBuilder(std::shared_ptr<Genotype> g, std::string name, std::vector<std::shared_ptr<Allele>>& alleles);
    explicit GenotypeBuilder(std::shared_ptr<Genotype> g);
    GenotypeBuilder& setAD(std::vector<int> AD);

    GenotypeBuilder& setDP(int DP);

    GenotypeBuilder& attributes(const std::map<std::string, AttributeValue>& attributes);

    GenotypeBuilder& attribute(const std::string& key, std::string& value);

    GenotypeBuilder& attribute(const std::string& key, int value);

    GenotypeBuilder& attribute(const std::string& key, std::vector<double>& value);

    GenotypeBuilder& attribute(const std::string& key, std::vector<int>& value);

    GenotypeBuilder& attribute(const std::string& key, std::vector<int>&& value);

    GenotypeBuilder buildPL(std::vector<double> GLs);
    std::shared_ptr<Genotype> make();

    /**
     * Set this genotype's alleles
     * @param alleles
     * @return
     */
    GenotypeBuilder& setAlleles(std::shared_ptr<std::vector<std::shared_ptr<Allele>>> alleles);

    GenotypeBuilder& setAlleles(std::vector<std::shared_ptr<Allele>>& alleles);

    GenotypeBuilder& phased(bool phased);
};


#endif //MUTECT2CPP_MASTER_GENOTYPEBUILDER_H
