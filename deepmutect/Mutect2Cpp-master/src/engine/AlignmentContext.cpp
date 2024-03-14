//
// Created by 梦想家xixi on 2022/1/11.
//

#include "AlignmentContext.h"

AlignmentContext::AlignmentContext(const std::list<pileRead *> &tumor, const std::list<pileRead *> &normal,
                                   int tumorSize, int normalSize, SimpleInterval &loc, int tid, SAMFileHeader *header)
		: tumor(tumor), normal(normal), loc(loc),
		  tid(tid), header(header), tumorSize(tumorSize), normalSize(normalSize) {
	// for debug
	// assert(tumorSize == tumor.size());
	// assert(normalSize == normal.size());
}

int AlignmentContext::getReadNum() const {
	return normalSize + tumorSize;
}

std::string AlignmentContext::getRefName() {
	return header->getSequenceDictionary().getSequences()[tid].getSequenceName();
}

int AlignmentContext::getPosition() const {
	return loc.getStart();
}

ReadPileup AlignmentContext::makeTumorPileup() {
	return {tid, loc.getStart(), tumor};
}

ReadPileup AlignmentContext::makeNormalPileup() {
	return {tid, loc.getStart(), normal};
}

bool AlignmentContext::isEmpty() const {
	return normalSize + tumorSize == 0;
}

SimpleInterval &AlignmentContext::getLocation() {
	return loc;
}
