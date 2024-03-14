//
// Created by 梦想家xixi on 2021/10/29.
//

#ifndef MUTECT2CPP_MASTER_ADAPTIVECHAINPRUNER_H
#define MUTECT2CPP_MASTER_ADAPTIVECHAINPRUNER_H

#include "Mutect2Utils.h"
#include <vector>
#include <set>
#include "./Path.h"
#include "ChainPruner.h"

template<class V, class E>
class AdaptiveChainPruner : public ChainPruner<V, E> {
private:
	double initialErrorProbability;
	double logOddsThreshold;
	int maxUnprunedVariants;

public:
	AdaptiveChainPruner(const double initialErrorProbability, const double logOddsThreshold,
	                    const int maxUnprunedVariants) : ChainPruner<V, E>(),
	                                                     initialErrorProbability(initialErrorProbability),
	                                                     logOddsThreshold(logOddsThreshold),
	                                                     maxUnprunedVariants(maxUnprunedVariants) {
		if (initialErrorProbability <= 0)
			throw std::invalid_argument("Must have positive error probability");
	}

protected:
	phmap::flat_hash_set<Path<V, E> *> chainsToRemove(std::vector<Path<V, E> *> chains) {
		if (chains.empty()) {
			phmap::flat_hash_set<Path<V, E> *> result;
			return result;
		}

		std::shared_ptr<DirectedSpecifics<V, E>> graph = chains[0]->getGraph();
		phmap::flat_hash_set<Path<V, E> *> probableErrorChains = likelyErrorChains(chains, graph, 0.001);
		int errorCount = 0, totalBases = 0;
		for (typename phmap::flat_hash_set<Path<V, E> *>::iterator siter = probableErrorChains.begin();
		     siter != probableErrorChains.end(); siter++) {
			errorCount += (*siter)->getLastEdge()->getMultiplicity();
		}
		for (typename std::vector<Path<V, E> *>::iterator viter = chains.begin(); viter != chains.end(); viter++) {
			for (typename std::vector<std::shared_ptr<E>>::iterator iter = (*viter)->getEdges().begin();
			     iter != (*viter)->getEdges().end(); iter++) {
				totalBases += (*iter)->getMultiplicity();
			}
		}
		double errorRate = (double) errorCount / totalBases;
		//std::cout << "probableErrorChains " << probableErrorChains.size() << std::endl;
		//std::cout << "errorCount " << errorCount << " totalBases " << totalBases << std::endl;
		return likelyErrorChains(chains, graph, errorRate);
	}

private:
	// Avoid getting bases multiple times
	struct chainWithBases {
		Path<V, E> *path = nullptr;
		std::shared_ptr<uint8_t[]> bases = nullptr;
		int len = 0;
		bool alreadyGotBases = false;

		explicit chainWithBases(Path<V, E> *path) : path(path) {}

		//.sorted(Comparator.comparingDouble((ToDoubleFunction<Path<V, E>>) chainLogOdds::get)
		//                        .reversed().thenComparingInt(Path::length))
		//according to JAVA version, chainLogOdds in descending order, if chainsLogOdds equal, length ascending order
		friend bool operator<(chainWithBases &a, chainWithBases &b) {
			if (a.path->getLogOdds() == b.path->getLogOdds()) {
				if (!a.alreadyGotBases) {
					a.bases = a.path->getBases(a.len);
					a.alreadyGotBases = true;
				}

				if (!b.alreadyGotBases) {
					b.bases = b.path->getBases(b.len);
					b.alreadyGotBases = true;
				}

				if (a.len != b.len)
					return a.len > b.len;

				return memcmp(a.bases.get(), b.bases.get(), a.len) < 0;
			}
			return a.path->getLogOdds() > b.path->getLogOdds();
		}
	};

	phmap::flat_hash_set<Path<V, E> *>
	likelyErrorChains(std::vector<Path<V, E> *> &chains, std::shared_ptr<DirectedSpecifics<V, E>> graph,
	                  double errorRate) {
		typename std::vector<Path<V, E> *>::iterator viter;
		phmap::flat_hash_set<Path<V, E> *> result;
		result.reserve(chains.size());

		for (viter = chains.begin(); viter != chains.end(); viter++) {
			(*viter)->setLogOdds(chainLogOdds(*viter, graph, errorRate));
			if ((*viter)->getLogOdds() < 2.302585092994046) {
				result.insert(*viter);
			}
		}

		std::vector<chainWithBases> newchains;
		for (viter = chains.begin(); viter != chains.end(); viter++) {
			if (isChainPossibleVariant(*viter, graph))
				newchains.template emplace_back(*viter);
		}
		std::sort(newchains.begin(), newchains.end());

		//maxUnprunedVariants = 100
		if (newchains.size() > 100) {
			for (auto iter = newchains.begin() + 100; iter != newchains.end(); iter++) {
				result.insert(iter->path);
			}
		}
		return result;
	}

	double chainLogOdds(Path<V, E> *chain, std::shared_ptr<DirectedSpecifics<V, E>> graph, double errorRate) {
		for (typename std::vector<std::shared_ptr<E>>::iterator viter = chain->getEdges().begin();
		     viter != chain->getEdges().end(); viter++) {
			if ((*viter)->getIsRef())
				return POSITIVE_INFINITY;
		}
		int leftTotalMultiplicity = 0, rightTotalMultiplicity = 0;
		phmap::flat_hash_set<std::shared_ptr<E>> outgoing = graph->outgoingEdgesOf(chain->getFirstVertex());
		phmap::flat_hash_set<std::shared_ptr<E>> incoming = graph->incomingEdgesOf(chain->getLastVertex());
		typename phmap::flat_hash_set<std::shared_ptr<E>>::iterator eiter;
		for (eiter = outgoing.begin(); eiter != outgoing.end(); eiter++) {
			leftTotalMultiplicity += (*eiter)->getMultiplicity();
		}
		for (eiter = incoming.begin(); eiter != incoming.end(); eiter++) {
			rightTotalMultiplicity += (*eiter)->getMultiplicity();
		}
		int leftMultiplicity = (chain->getEdges()[0])->getMultiplicity();
		int rightMultiplicity = chain->getLastEdge()->getMultiplicity();

		double leftLogOdds = graph->isSource(chain->getFirstVertex()) ? 0.0 : Mutect2Utils::logLikelihoodRatio(
				leftTotalMultiplicity - leftMultiplicity, leftMultiplicity, errorRate);
		double rightLogOdds = graph->isSink(chain->getLastVertex()) ? 0.0 : Mutect2Utils::logLikelihoodRatio(
				rightTotalMultiplicity - rightMultiplicity, rightMultiplicity, errorRate);

		return std::max(leftLogOdds, rightLogOdds);
	}

	bool isChainPossibleVariant(Path<V, E> *chain, std::shared_ptr<DirectedSpecifics<V, E>> &graph) {
		/*int len;
		std::shared_ptr<uint8_t[]> bases = chain->getBases(len);
		std::string baseString = std::string((char *) bases.get(), 0, len);
		if (baseString == "CCTTTACCCCTTTCAGCGATGTCCATTTTGTAA" || baseString == "TCTGTGAAGAGAAATGTACCCAGATCTATCATT") {
			std::cout << baseString <<" found\n\n";
		}*/
		typename phmap::flat_hash_set<std::shared_ptr<E>>::iterator eiter;
		int leftTotalMultiplicity = 0, rightTotalMultiplicity = 0;
		phmap::flat_hash_set<std::shared_ptr<E>> outgoing = graph->outgoingEdgesOf(chain->getFirstVertex());
		phmap::flat_hash_set<std::shared_ptr<E>> incoming = graph->incomingEdgesOf(chain->getLastVertex());
		for (eiter = outgoing.begin(); eiter != outgoing.end(); eiter++) {
			leftTotalMultiplicity += (*eiter)->getMultiplicity();
		}
		for (eiter = incoming.begin(); eiter != incoming.end(); eiter++) {
			rightTotalMultiplicity += (*eiter)->getMultiplicity();
		}
		int leftMultiplicity = (chain->getEdges()[0])->getMultiplicity();
		int rightMultiplicity = chain->getLastEdge()->getMultiplicity();

		return leftMultiplicity <= leftTotalMultiplicity / 2 || rightMultiplicity <= rightTotalMultiplicity / 2;
	}
};


#endif //MUTECT2CPP_MASTER_ADAPTIVECHAINPRUNER_H
