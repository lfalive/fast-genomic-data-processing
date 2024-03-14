//
// Created by 梦想家xixi on 2022/1/11.
//

#ifndef MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
#define MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H

#include "utils/ReadPileup.h"
#include "SimpleInterval.h"


class AlignmentContext {
private:
	const std::list<pileRead *> &tumor;
	const std::list<pileRead *> &normal;
	int tumorSize;
	int normalSize;
	SimpleInterval loc;
	int tid;
	SAMFileHeader *header;

public:
	AlignmentContext(const std::list<pileRead *> &tumor, const std::list<pileRead *> &normal, int tumorSize,
	                 int normalSize, SimpleInterval &loc, int tid, SAMFileHeader *header);

	AlignmentContext() = delete;

	int getReadNum() const;

	std::string getRefName();

	int getPosition() const;

	ReadPileup makeTumorPileup();

	ReadPileup makeNormalPileup();

	bool isEmpty() const;

	SimpleInterval &getLocation();
};


#endif //MUTECT2CPP_MASTER_ALIGNMENTCONTEXT_H
