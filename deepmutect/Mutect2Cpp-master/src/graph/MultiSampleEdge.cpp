//
// Created by 梦想家xixi on 2021/10/20.
//

#include "MultiSampleEdge.h"
#include <stdexcept>

MultiSampleEdge::MultiSampleEdge(bool isRef, int multiplicity, int singleSampleCapacity)
		: BaseEdge(isRef, multiplicity), singleSampleCapacity(singleSampleCapacity) {
	if (singleSampleCapacity <= 0) {
		throw std::invalid_argument("singleSampleCapacity must be > 0");
	}
	singleSampleMultiplicities.push(multiplicity);
	currentSingleSampleMultiplicity = multiplicity;
}

void MultiSampleEdge::flushSingleSampleMultiplicity() {
	singleSampleMultiplicities.push(currentSingleSampleMultiplicity);
	if (singleSampleMultiplicities.size() == singleSampleCapacity + 1) {
		singleSampleMultiplicities.pop();
	} else if (singleSampleMultiplicities.size() > singleSampleCapacity + 1) {
		throw std::invalid_argument("Somehow the per sample multiplicity list has grown too big");
	}
	currentSingleSampleMultiplicity = 0;
}

void MultiSampleEdge::incMultiplicity(int incr) {
	BaseEdge::incMultiplicity(incr);
	currentSingleSampleMultiplicity += incr;
}

int MultiSampleEdge::getPruningMultiplicity() const {
	return singleSampleMultiplicities.top();
}

bool MultiSampleEdge::operator==(const MultiSampleEdge &other) const {
	return (long) this == (long) &other;
}

bool MultiSampleEdge::operator<(const MultiSampleEdge &other) const {
	return (long) this < (long) &other;
}
