//
// Created by 梦想家xixi on 2021/11/30.
//

#include "GenotypeBuilder.h"
#include <utility>

GenotypeBuilder::GenotypeBuilder(std::shared_ptr<Genotype> g): sampleName(g->getSampleName()), alleles(g->getAlleles()), isPhased(g->isPhased()), GQ(g->getGQ()),  DP(g->getDP()), filters(g->getFilters())
{
    AD = g->getAD();
    PL = g->getPL();
    attributes(g->getExtendedAttributes());
}

GenotypeBuilder::GenotypeBuilder(std::shared_ptr<Genotype> g, std::string name, std::vector<std::shared_ptr<Allele>> &alleles):
    sampleName(name), alleles(alleles), isPhased(g->isPhased()), GQ(g->getGQ()), DP(g->getDP()), filters(g->getFilters())
{
    AD = g->getAD();
    PL = g->getPL();
    attributes(g->getExtendedAttributes());
}

std::shared_ptr<Genotype> GenotypeBuilder::create(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles) {
    return GenotypeBuilder(std::move(sampleName), std::move(alleles)).make();
}

std::shared_ptr<Genotype> GenotypeBuilder::make() {
    return std::make_shared<FastGenotype>(sampleName, alleles, isPhased, GQ, DP, AD, PL, filters, extendedAttributes);
}

std::shared_ptr<Genotype> GenotypeBuilder::create(std::string sampleName, std::vector<std::shared_ptr<Allele>> alleles,
                                  const std::map<std::string, AttributeValue>& attributes) {
    return GenotypeBuilder(std::move(sampleName), std::move(alleles)).attributes(attributes).make();
}

GenotypeBuilder& GenotypeBuilder::attributes(const std::map<std::string, AttributeValue>& attributes) {
    for(auto pairToAdd : attributes) {
        extendedAttributes.insert(pairToAdd);
    }
    return *this;
}

GenotypeBuilder& GenotypeBuilder::attribute(const std::string &key, std::string &value) {
    extendedAttributes.emplace(key, AttributeValue(value));
    return *this;
}

GenotypeBuilder& GenotypeBuilder::attribute(const std::string &key, int value) {
    extendedAttributes.emplace(key, AttributeValue(value));
    return *this;
}

GenotypeBuilder& GenotypeBuilder::attribute(const std::string &key, std::vector<double>& value){
    extendedAttributes.emplace(key, AttributeValue(value));
    return *this;
}

GenotypeBuilder& GenotypeBuilder::attribute(const std::string &key, std::vector<int>& value){
    extendedAttributes.emplace(key, AttributeValue(value));
    return *this;
}

GenotypeBuilder& GenotypeBuilder::attribute(const std::string &key, std::vector<int>&& value){
    extendedAttributes.emplace(key, AttributeValue(value));
    return *this;
}

std::shared_ptr<Genotype> GenotypeBuilder::create(const std::string& sampleName, const std::vector<std::shared_ptr<Allele>>& alleles, std::vector<double> gls) {
    return GenotypeBuilder(sampleName, alleles).buildPL(gls).make();
}

GenotypeBuilder GenotypeBuilder::buildPL(std::vector<double> GLs) {
    PL = GenotypeLikelihoods::fromLog10Likelihoods(GLs)->getAsPLs();
    return *this;
}

GenotypeBuilder&  GenotypeBuilder::setAD(std::vector<int> AD) {
    this->AD = std::move(AD);
    return *this;
}

GenotypeBuilder &GenotypeBuilder::setDP(int DP) {
    this->DP = DP;
    return *this;
}

GenotypeBuilder &GenotypeBuilder::setAlleles(std::shared_ptr<std::vector<std::shared_ptr<Allele>>> alleles) {
    if(alleles == nullptr)
        this->alleles.clear();
    else
        this->alleles = *alleles;
    return *this;
}

GenotypeBuilder& GenotypeBuilder::setAlleles(std::vector<std::shared_ptr<Allele>> &alleles) {
    if(alleles.empty())
        this->alleles.clear();
    else
        this->alleles = alleles;
    return *this;
}

GenotypeBuilder &GenotypeBuilder::phased(bool phased) {
    this->isPhased = phased;
    return *this;
}
