//
// Created by 梦想家xixi on 2021/11/27.
//

#ifndef MUTECT2CPP_MASTER_GENOTYPELIKELIHOODS_H
#define MUTECT2CPP_MASTER_GENOTYPELIKELIHOODS_H

#include <string>
#include <map>
#include <utility>
#include "GenotypeLikelihoodsAllelePair.h"
#include "Genotype.h"
#include "VariantContext.h"
#include "GenotypeType.h"

class Genotype;

class GenotypeLikelihoods {
private:
	static int numLikelihoodCache[5][10];
	std::vector<double> log10Likelihoods;
	std::string likelihoodsAsString_PLs;
	static GenotypeLikelihoodsAllelePair **diploidPLIndexToAlleleIndex;
	static int allelePairLength;

	GenotypeLikelihoods(std::string asString) : likelihoodsAsString_PLs(std::move(asString)),
	                                            log10Likelihoods(std::vector<double>()) {}

	GenotypeLikelihoods(std::vector<double> asVector) : log10Likelihoods(asVector) {}

	static std::vector<double> PLsToGLs(std::vector<int> pls);

	static std::vector<int> GLsToPLs(std::vector<double> GLs);

	static std::vector<double> parsePLsIntoLikelihoods(std::string likelihoodsAsString_PLs);

	static double maxPL(const std::vector<double> &GLs);

	static GenotypeLikelihoodsAllelePair **calculateDiploidPLcache(int altAlleles, int &length);

	static int calcNumLikelihoods(int FnumAlleles, int ploidy);

public:
	static GenotypeLikelihoods *fromLog10Likelihoods(std::vector<double> log10Likelihoods);

	static GenotypeLikelihoods *fromPLs(std::vector<int> pls);

	std::vector<double> getAsVector();

	std::vector<int> getAsPLs();

	static int calculatePLindex(int allele1Index, int allele2Index);

	static int numLikelihoods(int numAlleles, int ploidy);

	static GenotypeLikelihoodsAllelePair *getAllelePair(int PLindex);

	static std::vector<int> getAlleles(int PLindex, int ploidy);

	static void initial();

protected:
	static std::map<int, std::vector<std::vector<int>>> anyploidPloidyToPLIndexToAlleleIndices;
};


#endif //MUTECT2CPP_MASTER_GENOTYPELIKELIHOODS_H
