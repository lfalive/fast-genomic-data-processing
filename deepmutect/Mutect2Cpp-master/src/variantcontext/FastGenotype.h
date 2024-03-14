//
// Created by 梦想家xixi on 2021/11/30.
//

#ifndef MUTECT2CPP_MASTER_FASTGENOTYPE_H
#define MUTECT2CPP_MASTER_FASTGENOTYPE_H

#include "Genotype.h"
#include <iostream>

class FastGenotype : public Genotype {
private:
	std::vector<std::shared_ptr<Allele>> alleles;
	bool _isPhased;
	int GQ;
	int DP;
	std::vector<int> AD;
	std::vector<int> PL;
	std::map<std::string, AttributeValue> extendedAttributes;

public:
	FastGenotype(std::string sampleName, std::vector<std::shared_ptr<Allele>> &alleles, bool isPhased, int GQ, int DP,
	             std::vector<int> AD, std::vector<int> PL, const std::string &filters,
	             std::map<std::string, AttributeValue> extendedAttributes) :
			Genotype(std::move(sampleName), filters), alleles(alleles), _isPhased(isPhased), GQ(GQ), DP(DP),
			AD(std::move(AD)), PL(std::move(PL)), extendedAttributes(std::move(extendedAttributes)) {};

	std::vector<std::shared_ptr<Allele>> &getAlleles() override { return alleles; }

	std::shared_ptr<Allele> getAllele(int i) override { return alleles.at(i); }

	bool isPhased() override { return _isPhased; }

	int getDP() override { return DP; }

	std::vector<int> getAD() override { return AD; };

	int getGQ() override { return GQ; }

	std::vector<int> getPL() override { return PL; };

	std::map<std::string, AttributeValue> &getExtendedAttributes() { return extendedAttributes; }
};


#endif //MUTECT2CPP_MASTER_FASTGENOTYPE_H
