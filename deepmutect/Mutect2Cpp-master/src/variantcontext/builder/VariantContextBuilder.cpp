//
// Created by 梦想家xixi on 2021/12/6.
//

#include "VariantContextBuilder.h"

#include <memory>
#include <utility>

VariantContextBuilder::VariantContextBuilder(std::string & source, std::string & contig, long start, long stop,
                                             std::shared_ptr<std::vector<std::shared_ptr<Allele>>>   alleles) : genotypes(nullptr), log10PError(1.0),
                                             attributesCanBeModified(false), source(source), contig(contig), start(start), stop(stop), alleles(std::move(alleles)), filters(
                nullptr), attribute(nullptr){
    toValidate.insert(ALLELES);
}

std::shared_ptr<VariantContext> VariantContextBuilder::make(bool leaveModifyableAsIs) {
    if(!leaveModifyableAsIs) {
        attributesCanBeModified = false;
    }
    return std::make_shared<VariantContext>(source, ID, contig, start, stop, alleles, genotypes, log10PError, filters, attribute, fullyDecoded,toValidate);
}

std::shared_ptr<VariantContext> VariantContextBuilder::make() {
    return make(false);
}

VariantContextBuilder::VariantContextBuilder(std::shared_ptr<VariantContext> &parent) : alleles(std::make_shared<std::vector<std::shared_ptr<Allele>>>(parent->getAlleles()) ), attribute(parent->getAttributesAsPointer()), attributesCanBeModified(false),
contig(parent->getContig()), filters(parent->getFiltersMaybeNull()), genotypes(parent->getGenotypes()), ID(parent->getID()), log10PError(parent->getLog10PError()),
source(parent->getSource()), start(parent->getStart()), stop(parent->getEnd()), fullyDecoded(parent->isFullyDecoded()){}

void VariantContextBuilder::setStart(long start) {
    this->start = start;
}

void VariantContextBuilder::setStop(long stop) {
    this->stop = stop;
}

VariantContextBuilder* VariantContextBuilder::setAlleles(const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> &  alleles) {
    this->alleles = alleles;
    toValidate.insert(ALLELES);
    return this;
}

VariantContextBuilder* VariantContextBuilder::setAlleles(phmap::flat_hash_set<std::shared_ptr<Allele>> &alleles) {
    this->alleles = std::make_shared<std::vector<std::shared_ptr<Allele>>>();
    for(const auto & allele : alleles)
    {
        this->alleles->emplace_back(allele);
    }
    toValidate.insert(ALLELES);
    return this;
}

VariantContextBuilder::~VariantContextBuilder() {
    alleles->clear();
}

VariantContextBuilder::VariantContextBuilder(): genotypes(nullptr), filters(nullptr), attribute(nullptr), log10PError(0.0) {

}

VariantContextBuilder *VariantContextBuilder::setSource(std::string& source) {
    this->source = source;
    return this;
}

VariantContextBuilder *VariantContextBuilder::setId(std::string &ID) {
    this->ID = ID;
    return this;
}

VariantContextBuilder *VariantContextBuilder::setLoc(std::string &contig, long start, long stop) {
    this->contig = contig;
    this->start = start;
    this->stop = stop;
    toValidate.insert(ALLELES);
    return this;
}

VariantContextBuilder *VariantContextBuilder::setGenotypes(std::shared_ptr<GenoTypesContext> genotypesContext) {
    this->genotypes = genotypesContext;
    if(genotypes != nullptr)
    {
        genotypes->setImmutable();
        toValidate.insert(GENOTYPES);
    }
    return this;
}

VariantContextBuilder *VariantContextBuilder::setGenotypes(std::vector<std::shared_ptr<Genotype>>& genotypes) {
    return setGenotypes(std::make_shared<GenoTypesContext>(genotypes));
}

VariantContextBuilder *VariantContextBuilder::setLog10PError(double log10PError) {
    this->log10PError = log10PError;
    return this;
}

VariantContextBuilder *VariantContextBuilder::setFilters(std::set<std::string> *filters) {
    if(!filters)
    {
        unfiltered();
    } else {
        this->filtersCanBeModified = false;
        this->filters = filters;
        toValidate.insert(FILTERS);
    }
    return this;
}

VariantContextBuilder *VariantContextBuilder::unfiltered() {
    this->filters = nullptr;
    this->filtersCanBeModified = false;
    return this;
}

VariantContextBuilder *
VariantContextBuilder::setAttributes(std::shared_ptr<std::map<std::string, AttributeValue>> attribute) {
    this->attribute = attribute;
    this->attributesCanBeModified = true;
    return this;
}

VariantContextBuilder *VariantContextBuilder::setAttributes(std::shared_ptr<std::map<std::string, std::vector<double>>> attribute) {
    this->attribute = std::make_shared<std::map<std::string, AttributeValue>>();
    if(attribute != nullptr)
    {
        for(auto & iter : *attribute)
        {
            this->attribute->emplace(iter.first, AttributeValue(iter.second));
        }
    }
    this->attributesCanBeModified = true;
    return this;
}

VariantContextBuilder *VariantContextBuilder::setAttribute(const std::string &key, int value) {
    this->attribute->emplace(key, AttributeValue(value));
    return this;
}

VariantContextBuilder *VariantContextBuilder::setAttribute(const std::string &key, std::vector<double> value) {
    this->attribute->emplace(key, AttributeValue(value));
    return this;
}

VariantContextBuilder* VariantContextBuilder::setAttribute(const std::string& key, std::shared_ptr<std::vector<double>> value)
{
    return setAttribute(key, *value);
}