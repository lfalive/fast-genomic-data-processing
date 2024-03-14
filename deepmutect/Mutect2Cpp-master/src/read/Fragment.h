//
// Created by lhh on 5/20/22.
//

#ifndef MUTECT2CPP_MASTER_FRAGMENT_H
#define MUTECT2CPP_MASTER_FRAGMENT_H


#include <Locatable.h>
#include "SimpleInterval.h"
#include "samtools/SAMRecord.h"

class Fragment : public Locatable{
private:
    SimpleInterval interval;
    std::vector<std::shared_ptr<SAMRecord>> reads;

public:
    Fragment(std::shared_ptr<SAMRecord>);

    Fragment(std::vector<std::shared_ptr<SAMRecord>>& reads);

    static std::shared_ptr<Fragment> create(std::vector<std::shared_ptr<SAMRecord>>& reads);

    static std::shared_ptr<Fragment> create(std::vector<std::shared_ptr<SAMRecord>>&& reads);

    static std::shared_ptr<Fragment> createAndAvoidFailure(std::vector<std::shared_ptr<SAMRecord>>& reads);

    std::string getContig() const;

    int getStart() const;
    int getEnd() const;
};


#endif //MUTECT2CPP_MASTER_FRAGMENT_H
