//
// Created by lhh on 11/12/21.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPE_H
#define MUTECT2CPP_MASTER_GENOTYPE_H

#include "vector"
#include "string"
#include "Allele.h"
#include "GenotypeLikelihoods.h"
#include "GenotypeType.h"
#include "AttributeValue.h"


class GenotypeLikelihoods;

class Genotype {
public:
    static const long serialVersionUID = 1L;
    const static std::string PRIMARY_KEYS[6];
    const static std::string PHASED_ALLELE_SEPARATOR;
    const static std::string UNPHASED_ALLELE_SEPARATOR;

    virtual std::vector<std::shared_ptr<Allele>> &getAlleles() = 0;

    virtual std::shared_ptr<Allele> getAllele(int var1) = 0;

    virtual bool isPhased() = 0;

    int getPloidy() { return getAlleles().size(); }

    virtual int getDP() = 0;

    virtual std::vector<int> getAD() = 0;

    virtual int getGQ() = 0;

    virtual std::vector<int> getPL() = 0;

    std::string getSampleName() const { return sampleName; }

    int countAllele(const std::shared_ptr<Allele> &allele);

    bool hasPL();

    bool hasAD();

    bool hasGQ() { return this->getGQ() != -1; }

    bool hasDP() { return this->getDP() != -1; }

    bool isHom() { return isHomRef() || isHomVar(); }

    bool isHomRef() { return getType() == HOM_REF; }

    bool isHomVar() { return getType() == HOM_VAR; }

    GenotypeType getType();

    bool isHet() { return getType() == HET; }

    bool isHetNonRef();

    bool isNoCall() { return getType() == NO_CALL; }

    bool isCalled() { return getType() != NO_CALL && getType() != UNAVAILABLE; }

    bool isMixed() { return getType() == MIXED; }

    bool isAvailable() { return getType() != UNAVAILABLE; }

    bool hasLikelihoods();

    GenotypeLikelihoods *getLikelihoods();

    bool isNonInformative();

    std::string getGenotypeString(bool ignoreRefState);

    bool operator<(const Genotype &other) const;

    bool sameGenotype(Genotype *other, bool ignorePhase);

    bool sameGenotype(Genotype *other);

    virtual std::map<std::string, AttributeValue> &getExtendedAttributes() = 0;

    bool hasExtendedAttribute(const std::string &key);

    AttributeValue getExtendedAttribute(const std::string &, void *defaultValue);

    AttributeValue getExtendedAttribute(const std::string &);

    std::string getFilters() { return filters; }

    bool isFiltered() { return !filters.empty(); }

    bool hasAnyAttribute(std::string key);

private:
    std::string sampleName;
    GenotypeType type = GenotypeType_NULL;
    std::string filters;

protected:
    Genotype(std::string sampleName, const std::string &filters);

    GenotypeType determineType();

    std::vector<std::string> getAlleleStrings();
};


#endif //MUTECT2CPP_MASTER_GENOTYPE_H
