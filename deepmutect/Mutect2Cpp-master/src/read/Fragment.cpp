//
// Created by lhh on 5/20/22.
//

#include <cassert>
#include "Fragment.h"

Fragment::Fragment(std::shared_ptr<SAMRecord> read): interval(read->getAssignedContig(), std::min(read->getStart(), read->getEnd()),
                                                              std::max(read->getStart(), read->getEnd()))
{
    reads.push_back(read);
}

Fragment::Fragment(std::vector<std::shared_ptr<SAMRecord>>& reads)
{
    assert(reads.size() == 2);
    this->reads = reads;
    int start = std::min(reads[0]->getStart(), reads[1]->getStart());
    int end = std::max(reads[0]->getEnd(), reads[1]->getEnd());
    interval = SimpleInterval(reads[0]->getContig(), std::min(start, end), std::max(start, end));
}

std::shared_ptr<Fragment> Fragment::create(std::vector<std::shared_ptr<SAMRecord>> &reads) {
    assert(reads.size() <= 2);
    assert(!reads.empty());
    return reads.size() == 1 ? std::make_shared<Fragment>(reads[0]) : std::make_shared<Fragment>(reads);
}

std::shared_ptr<Fragment> Fragment::create(std::vector<std::shared_ptr<SAMRecord>> &&reads) {
    assert(reads.size() <= 2);
    assert(!reads.empty());
    return reads.size() == 1 ? std::make_shared<Fragment>(reads[0]) : std::make_shared<Fragment>(reads);
}

std::shared_ptr<Fragment> Fragment::createAndAvoidFailure(std::vector<std::shared_ptr<SAMRecord>>& reads) {
    if(reads.size() <= 2)
        return create(reads);
    else {
        std::vector<std::shared_ptr<SAMRecord>> nonSupplementaryReads;
        for(auto & read : reads)
        {
            if(!(read->isDuplicate() || read->isSecondaryAlignment() || read->isSupplementaryAlignment()))
                nonSupplementaryReads.emplace_back(read);
        }
        if(nonSupplementaryReads.size() > 2)
        {
            std::cerr << "Fragment::createAndAvoidFailure() " << "More than two reads with the same name found.  Using two reads randomly to combine as a fragment.";
            std::vector<std::shared_ptr<SAMRecord>> temp{reads[0], reads[1]};
            return create(temp);
        }
        if(nonSupplementaryReads.empty())
            return create({reads[0]});
        return create(nonSupplementaryReads);
    }
}

std::string Fragment::getContig() const
{
    return interval.getContig();
}

int Fragment::getStart() const
{
    return interval.getStart();
}

int Fragment::getEnd() const
{
    return interval.getEnd();
}