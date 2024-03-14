//
// Created by 梦想家xixi on 2021/11/27.
//

#include "GenotypeLikelihoods.h"
#include <regex>
#include <vector>
#include <cmath>
#include <sstream>
#include "StringUtils.h"

int GenotypeLikelihoods::numLikelihoodCache[5][10] = {0};

int GenotypeLikelihoods::allelePairLength = 0;

GenotypeLikelihoodsAllelePair **GenotypeLikelihoods::diploidPLIndexToAlleleIndex = nullptr;

std::map<int, std::vector<std::vector<int>>> GenotypeLikelihoods::anyploidPloidyToPLIndexToAlleleIndices;

GenotypeLikelihoods *GenotypeLikelihoods::fromLog10Likelihoods(std::vector<double> log10Likelihoods) {
	return new GenotypeLikelihoods(log10Likelihoods);
}

std::vector<double> GenotypeLikelihoods::PLsToGLs(std::vector<int> pls) {
	std::vector<double> likelihoodsAsVector(pls.size());
	for (int i = 0; i < pls.size(); i++) {
		likelihoodsAsVector[i] = (double) pls[i] / -10.0;
	}
	return likelihoodsAsVector;
}

GenotypeLikelihoods *GenotypeLikelihoods::fromPLs(std::vector<int> pls) {
	return new GenotypeLikelihoods(PLsToGLs(pls));
}

std::vector<double> GenotypeLikelihoods::getAsVector() {
	if (log10Likelihoods.empty()) {
		log10Likelihoods = parsePLsIntoLikelihoods(likelihoodsAsString_PLs);
	}
	return log10Likelihoods;
}

std::vector<double> GenotypeLikelihoods::parsePLsIntoLikelihoods(std::string likelihoodsAsString_PLs) {
	if (likelihoodsAsString_PLs != ".") {
		std::regex ws_re(",");
		std::vector<std::string> strings(
				std::sregex_token_iterator(likelihoodsAsString_PLs.begin(), likelihoodsAsString_PLs.end(), ws_re, -1),
				std::sregex_token_iterator());
		std::vector<double> likelihoodsAsVector(strings.size());
		for (int i = 0; i < strings.size(); i++) {
			likelihoodsAsVector[i] = StringUtils::parseInt(strings[i]) / -10.0;
		}
		return likelihoodsAsVector;
	} else
		return {};
}

std::vector<int> GenotypeLikelihoods::getAsPLs() {
	std::vector<double> GLs = getAsVector();
	return GLs.empty() ? std::vector<int>() : GLsToPLs(GLs);
}

std::vector<int> GenotypeLikelihoods::GLsToPLs(std::vector<double> GLs) {
	std::vector<int> pls(GLs.size());
	double adjust = maxPL(GLs);
	for (int i = 0; i < GLs.size(); ++i) {
		pls[i] = (int) std::round(std::min(-10.0 * (GLs[i] - adjust), 2.147483647E9));
	}
	return pls;
}

double GenotypeLikelihoods::maxPL(const std::vector<double> &GLs) {
	return std::max(-1.0 / 0.0, *max_element(GLs.begin(), GLs.end()));
}

int GenotypeLikelihoods::calculatePLindex(int allele1Index, int allele2Index) {
	return allele2Index * (allele2Index + 1) / 2 + allele1Index;
}

GenotypeLikelihoodsAllelePair **GenotypeLikelihoods::calculateDiploidPLcache(int altAlleles, int &length) {
	int numLikelihood = numLikelihoods(1 + altAlleles, 2);
	int i;
	GenotypeLikelihoodsAllelePair **cache = new GenotypeLikelihoodsAllelePair *[numLikelihood];
	for (i = 0; i <= altAlleles; ++i) {
		for (int allele2 = i; allele2 <= altAlleles; ++allele2) {
			cache[calculatePLindex(i, allele2)] = new GenotypeLikelihoodsAllelePair(i, allele2);
		}
	}

	for (i = 0; i < numLikelihood; ++i) {
		if (cache[i] == nullptr) {
			for (int k = 0; k < numLikelihood; ++k) {
				delete cache[k];
			}
			delete[] cache;
			throw std::invalid_argument("BUG: cache entry is unexpected null");
		}
	}
	length = numLikelihood;
	return cache;
}

int GenotypeLikelihoods::numLikelihoods(int numAlleles, int ploidy) {
	return numAlleles < 5 && ploidy < 10 ? numLikelihoodCache[numAlleles][ploidy]
	                                     : calcNumLikelihoods(numAlleles, ploidy);
}

int GenotypeLikelihoods::calcNumLikelihoods(int numAlleles, int ploidy) {
	if (numAlleles == 1) {
		return 1;
	} else if (ploidy == 1) {
		return numAlleles;
	} else {
		int acc = 0;
		for (int k = 0; k <= ploidy; ++k) {
			acc += calcNumLikelihoods(numAlleles - 1, ploidy - k);
		}
		return acc;
	}
}

GenotypeLikelihoodsAllelePair *GenotypeLikelihoods::getAllelePair(int PLindex) {
	if (PLindex >= 0 && PLindex < allelePairLength) {
		return diploidPLIndexToAlleleIndex[PLindex];
	} else {
		throw std::invalid_argument("The PL index cannot be negative");
	}
}

//TODO:需要同步
std::vector<int> GenotypeLikelihoods::getAlleles(int PLindex, int ploidy) {
	if (ploidy == 2) {
		GenotypeLikelihoodsAllelePair *tmpPair = getAllelePair(PLindex);
		return std::vector<int>{tmpPair->alleleIndex1, tmpPair->alleleIndex2};
	} else if (anyploidPloidyToPLIndexToAlleleIndices.find(ploidy) == anyploidPloidyToPLIndexToAlleleIndices.end()) {
		throw std::invalid_argument("Must initialize the cache of allele anyploid indices for ploidy ");
	} else if (PLindex >= 0 && PLindex < anyploidPloidyToPLIndexToAlleleIndices.at(ploidy).size()) {
		return anyploidPloidyToPLIndexToAlleleIndices.at(ploidy).at(PLindex);
	} else {
		throw std::invalid_argument("cannot have a negative value.");
	}
}

void GenotypeLikelihoods::initial() {
	for (int numAlleles = 1; numAlleles < 5; ++numAlleles) {
		for (int ploidy = 1; ploidy < 10; ++ploidy) {
			numLikelihoodCache[numAlleles][ploidy] = calcNumLikelihoods(numAlleles, ploidy);
		}
	}
	diploidPLIndexToAlleleIndex = calculateDiploidPLcache(50, allelePairLength);
}
