//
// Created by 梦想家xixi on 2021/10/20.
//

#include <stdexcept>
#include "BaseEdge.h"
#include "Mutect2Utils.h"
#include "graph/utils/GraphObjectPool.h"

BaseEdge::BaseEdge(bool isRef, int multiplicity) : isRef(isRef), multiplicity(multiplicity) {
	if (multiplicity < 0)
		throw std::invalid_argument("multiplicity must be >= 0");
}

void BaseEdge::incMultiplicity(const int incr) {
	if (incr < 0)
		throw std::invalid_argument("incr must be >= 0");
	multiplicity += incr;
}

void BaseEdge::setMultiplicity(int value) {
	if (value < 0)
		throw std::invalid_argument("multiplicity must be >= 0");
	multiplicity = value;
}

void BaseEdge::setIsRef(const bool ref) {
	this->isRef = ref;
}

BaseEdge BaseEdge::add(BaseEdge &edge) {
	multiplicity += edge.getMultiplicity();
	isRef = isRef || edge.getIsRef();
	return *this;
}

std::shared_ptr<BaseEdge> BaseEdge::makeOREdge(const std::vector<std::shared_ptr<BaseEdge>> &edges, int multiplicity) {
	if (edges.empty())
		throw std::invalid_argument("have no edge");
	bool anyRef = false;
	for (const std::shared_ptr<BaseEdge> &edge: edges) {
		if (edge->getIsRef()) {
			anyRef = true;
			break;
		}
	}
	return GraphObjectPool::createSeqEdge(anyRef, multiplicity);
}

std::shared_ptr<BaseEdge>
BaseEdge::makeOREdge(const phmap::flat_hash_set<std::shared_ptr<BaseEdge>> &edges, int multiplicity) {
	if (edges.empty())
		throw std::invalid_argument("have no edge");
	bool anyRef = false;
	for (const std::shared_ptr<BaseEdge> &edge: edges) {
		if (edge->getIsRef()) {
			anyRef = true;
			break;
		}
	}
	return GraphObjectPool::createSeqEdge(anyRef, multiplicity);
}

