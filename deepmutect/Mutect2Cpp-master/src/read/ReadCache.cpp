//
// Created by 梦想家xixi on 2022/1/11.
//

#include <cassert>
#include <random>
#include "ReadCache.h"
#include "ReadUtils.h"

#define START_GAP 50

// unused constructor
ReadCache::ReadCache(aux_t **data, std::vector<char*> & bam_name, std::shared_ptr<ReferenceCache> & cache) : data(data), tid(0), bam_name(bam_name),
                                                                    readTransformer(cache, data[0]->header, 5){
    bam1_t * b;
    b = bam_init1();
    std::string region(sam_hdr_tid2name(data[0]->hdr, 0));
    for(int i = 0; i < bam_name.size(); i++){
        int result;
        //int count = 0;
        hts_idx_t * idx = sam_index_load(data[i]->fp, bam_name[i]);
        if(idx == nullptr)
            throw std::invalid_argument("random alignment retrieval only works for indexed BAM or CRAM files.");
        hts_idxes.push_back(idx);
        hts_itr_t* iter = sam_itr_querys(idx, data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
            bam1_t * transformed_read = readTransformer.apply(b, data[i]->hdr);
            if(ReadFilter::test(transformed_read, data[i]->hdr)) {
                std::shared_ptr<SAMRecord> read = std::make_shared<SAMRecord>(transformed_read, data[i]->hdr);
                if(i == 0) {
                    read->setGroup(0);
                    normalReads.emplace(getpileRead(read));
                }
                else {
                    read->setGroup(1);
                    tumorReads.emplace(getpileRead(read));
                }
            }
            //if(count > 500)
            //    break;
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(b);

    //start = end = currentPose = 0;
}

ReadCache::ReadCache(aux_t **data, std::vector<char *> &bam_name, int tid, int start, int end, int maxReadsPerAlignmentStart, int maxAssemblyRegionSize, std::shared_ptr<ReferenceCache> & cache, bool bqsr_within_mutect, BQSRReadTransformer * tumorTransformer, BQSRReadTransformer * normalTransformer) :
tid(tid), data(data), bam_name(bam_name), readTransformer(cache, data[0]->header, 5), bqsr_within_mutect(bqsr_within_mutect), tumorTransformer(tumorTransformer), normalTransformer(normalTransformer){

    for(int i = 0; i < bam_name.size(); i++){
        hts_idx_t * idx = sam_index_load(data[i]->fp, bam_name[i]);
        if(idx == nullptr)
        {
            throw std::invalid_argument("random alignment retrieval only works for indexed BAM or CRAM files.");
            exit(1);
        }

        hts_idxes.push_back(idx);
    }

    chr_len = sam_hdr_tid2len(data[0]->hdr, tid);
    ExtendedStart = max(start - maxAssemblyRegionSize, 0);
    ExtendedEnd = min(end + maxAssemblyRegionSize, chr_len);
    chr_name = std::string(sam_hdr_tid2name(data[0]->hdr, tid));
	chr_name_int = ContigMap::getContigInt(chr_name);
    currentPose = ExtendedStart - 1;

	if (maxReadsPerAlignmentStart < 0)
		maxCoverage = DEFAULT_MAX_COVERAGE;
	else if (maxReadsPerAlignmentStart == 0)
		maxCoverage = INT_MAX;
	else
		maxCoverage = maxReadsPerAlignmentStart;

    readData(string(sam_hdr_tid2name(data[0]->hdr, tid)) + ":" + to_string(ExtendedStart) + "-" + to_string(ExtendedEnd));

    if(bqsr_within_mutect)
    {
        assert(tumorTransformer);
        assert(normalTransformer);
    }
}

ReadCache::~ReadCache() {
    clear();

    for(auto idx : hts_idxes)
    {
        hts_idx_destroy(idx);
    }
}

void ReadCache::readData(const string &region)
{
    bam1_t * b;
    b = bam_init1();
	std::vector<pileRead *> pendingReads;

	for(int i = 0; i < bam_name.size(); i++){
        int result;
        hts_itr_t* iter = sam_itr_querys(hts_idxes[i], data[i]->hdr, region.c_str());
        while((result = sam_itr_next(data[i]->fp, iter, b)) >= 0) {
	        if(ReadFilter::test(b, data[i]->hdr)) {
				// recalibrate base qualities
				if(bqsr_within_mutect) {
					if(i == 0) {
						normalTransformer->apply(b);
					} else {
						tumorTransformer->apply(b);
					}
				}

				bam1_t * transformed_read = readTransformer.apply(b, data[i]->hdr);
				std::shared_ptr<SAMRecord> read = std::make_shared<SAMRecord>(transformed_read, data[i]->hdr);
				read->setGroup(i);
				pendingReads.emplace_back(getpileRead(read));

				if(transformed_read != b)
					bam_destroy1(transformed_read);
			}
		}
        hts_itr_destroy(iter);
    }

	// downSample
	if (BOOST_LIKELY(!pendingReads.empty())) {
		// sort reads by coordinate
		std::sort(pendingReads.begin(), pendingReads.end(), [](pileRead *a, pileRead *b) -> bool {
			return a->read->getStart() < b->read->getStart();
		});
		int last_start = pendingReads.front()->read->getStart(), count = 1;
		for (int i = 1; i <= pendingReads.size(); ++i) {
			if (BOOST_LIKELY(i != pendingReads.size()) && pendingReads[i]->read->getStart() == last_start) {
				count++;
				continue;
			}

			// handle position change
			if (count > maxCoverage) {
				std::vector<pileRead *> toDownSample(pendingReads.begin() + i - count, pendingReads.begin() + i);
				std::vector<pileRead *> downSampled = downSample(toDownSample);
				// std::cout << (*toDownSample.begin())->read->getStart() + 1 << " " << toDownSample.size() << " " << downSampled.size() << std::endl;
				splitPendingReads(downSampled, 0 ,(int)downSampled.size());
			} else {
				splitPendingReads(pendingReads, i - count, i);
			}

			if (BOOST_UNLIKELY(i == pendingReads.size())) break;
			last_start = pendingReads[i]->read->getStart(), count = 1;
		}
	}

    bam_destroy1(b);
}

void ReadCache::splitPendingReads(const std::vector<pileRead *>& pendingReads, int start, int end) {
	for (int j = start; j < end; ++j) {
		pileRead* finalizedPileRead = pendingReads[j];
		if (finalizedPileRead->read->getGroup() == 0) {
			normalReads.emplace(finalizedPileRead);
		} else {
			tumorReads.emplace(finalizedPileRead);
		}
	}
}

std::vector<pileRead *> ReadCache::downSample(const std::vector<pileRead *>& pendingReads) const {
	std::vector<pileRead *> filtered;
	filtered.reserve(pendingReads.size());

	for (const auto &pileRead: pendingReads) {
		if (pileRead->read->getMappingQuality() > SUSPICIOUS_MAPPING_QUALITY) {
			filtered.emplace_back(pileRead);
		}
	}

	if (filtered.size() <= maxCoverage)
		return filtered;

	shuffle(filtered.begin(), filtered.end(), std::mt19937(std::random_device()()));
	return {filtered.begin(), filtered.begin() + maxCoverage};
}

int ReadCache::getNextPos() {
    int nextPose = currentPose + 1;
    if(nextPose > chr_len || nextPose > ExtendedEnd) {
        throw std::invalid_argument("please check first");
    }
    currentPose = nextPose;
    return nextPose;
}

bool ReadCache::hasNextPos() {
    int nextPose = currentPose + 1;
    return nextPose <= ExtendedEnd;
}

/*void ReadCache::advanceLoad() {
    // clear the elements of last region
    if(!normalReads.empty() && !tumorReads.empty())
        throw std::logic_error("error");
    clear();

    start = end + 1;
    end = start + REGION_SIZE -1;

    std::string region = std::string(sam_hdr_tid2name(data[0]->hdr, tid)) + ':' +
            std::to_string(start+1) + '-' + std::to_string(end);

    readData(region);
}*/

AlignmentContext ReadCache::getAlignmentContext() {
    getNextPos();
    // remove some read if necessary
    auto iter = tumorReadsForAlignment.begin();
    while(iter != tumorReadsForAlignment.end() && (*iter)->activateStop < currentPose) {
        delete *iter;
        tumorReadsForAlignment.erase(iter++);
		tumorReadsForAlignmentSize--;
    }
    iter = normalReadsForAlignment.begin();
    while(iter != normalReadsForAlignment.end() && (*iter)->activateStop < currentPose) {
        delete *iter;
        normalReadsForAlignment.erase(iter++);
		normalReadsForAlignmentSize--;
    }

    while(!tumorReads.empty() && (*tumorReads.front()).read->getStart() <= currentPose){
        if((*tumorReads.front()).activateStart <= currentPose) {
            if((*tumorReads.front()).activateStop >= currentPose) {
                InsertPileToAlignment(tumorReads.front(), tumorReadsForAlignment);
				tumorReadsForAlignmentSize++;
				tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                tumorReads.pop();
            } else {
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                delete tumorReads.front();
                tumorReads.pop();
            }
        } else {
            if(InsertPileToCache(tumorReads.front(), tumorCache)) {
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                tumorReads.pop();
            } else {
                tumorReadsForRegion.emplace_back(tumorReads.front()->read);
                delete tumorReads.front();
                tumorReads.pop();
            }
        }
    }
    while(!normalReads.empty() && (*normalReads.front()).read->getStart() <= currentPose){
        if((*normalReads.front()).activateStart <= currentPose) {
            if((*normalReads.front()).activateStop >= currentPose) {
                InsertPileToAlignment(normalReads.front(), normalReadsForAlignment);
	            normalReadsForAlignmentSize++;
				normalReadsForRegion.emplace_back(normalReads.front()->read);
                normalReads.pop();
            } else {
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                delete normalReads.front();
                normalReads.pop();
            }
        } else {
            if(InsertPileToCache(normalReads.front(), normalCache)) {
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                normalReads.pop();
            } else {
                normalReadsForRegion.emplace_back(normalReads.front()->read);
                delete normalReads.front();
                normalReads.pop();
            }
        }
    }
    iter = tumorCache.begin();
    while (iter != tumorCache.end() && (*iter)->activateStart <= currentPose) {
        InsertPileToAlignment(*iter, tumorReadsForAlignment);
	    tumorReadsForAlignmentSize++;
		tumorCache.erase(iter++);
    }
    iter = normalCache.begin();
    while (iter != normalCache.end() && (*iter)->activateStart <= currentPose) {
        InsertPileToAlignment(*iter, normalReadsForAlignment);
	    normalReadsForAlignmentSize++;
		normalCache.erase(iter++);
    }
    SimpleInterval loc(chr_name_int, currentPose, currentPose);
    return {tumorReadsForAlignment, normalReadsForAlignment, tumorReadsForAlignmentSize, normalReadsForAlignmentSize, loc, tid, data[0]->header};
}

std::vector<std::shared_ptr<SAMRecord>> ReadCache::getReadsForRegion(AssemblyRegion & region) {
    std::vector<std::shared_ptr<SAMRecord>> ret;
    std::shared_ptr<SimpleInterval> loc = region.getExtendedSpan();
    auto iter = tumorReadsForRegion.begin();
    while(iter != tumorReadsForRegion.end()) {
        std::shared_ptr<SimpleInterval> readLoc = (*iter)->getLoc();
        if((*iter)->getEndAfterFliter() < loc->getStart()) {
            tumorReadsForRegion.erase(iter++);
        } else {
            if(loc->overlaps(readLoc)) {
                ret.emplace_back((*iter));
            }
            iter++;
        }
    }

    iter = normalReadsForRegion.begin();
    while(iter != normalReadsForRegion.end()) {
        std::shared_ptr<SimpleInterval> readLoc = (*iter)->getLoc();
        if((*iter)->getEndAfterFliter() < loc->getStart()) {
            normalReadsForRegion.erase(iter++);
        } else {
            if(loc->overlaps(readLoc)) {
                ret.emplace_back((*iter));
            }
            iter++;
        }
    }
    return ret;
}

pileRead *ReadCache::getpileRead(const std::shared_ptr<SAMRecord> &read) {
    int adaptorBoundary = read->getAdaptorBoundary();
    if (adaptorBoundary == INT32_MIN || read->getFragmentLength() > 100)
        return new pileRead{read, read->getStart(), read->getEnd()};
    int start = 0;
    int end = 0;
    if(read->isReverseStrand()) {
        start = std::max(adaptorBoundary+1, read->getStart());
        end = read->getEnd();

    } else {
        start = read->getStart();
        end = std::min(adaptorBoundary-1, read->getEnd());
    }
    return new pileRead{read, start, end};
}

void ReadCache::InsertPileToAlignment(pileRead* pile, std::list<pileRead*> & toAdd) {
    auto iter = toAdd.end();
    iter--;
    while(iter != toAdd.begin())
    {
        if((*iter)->activateStop < pile->activateStop)
            break;
        iter--;
    }
    if(!toAdd.empty() && (*iter)->activateStop < pile->activateStop)
        toAdd.insert(++iter, pile);
    else {
        toAdd.emplace_front(pile);
    }
}

bool ReadCache::InsertPileToCache(pileRead *pile, std::list<pileRead *> & toAdd) {
    if(pile->activateStart > pile->activateStop)
        return false;
    auto iter = toAdd.end();
    iter--;
    while(iter != toAdd.begin())
    {
        if((*iter)->activateStart < pile->activateStart)
            break;
        iter--;
    }
    if(!toAdd.empty() && (*iter)->activateStart < pile->activateStart)
        toAdd.insert(++iter, pile);
    else {
        toAdd.emplace_front(pile);
    }
    return true;
}

void ReadCache::clear()
{
    normalReadsForRegion.clear();
    tumorReadsForRegion.clear();

    for(pileRead* read : tumorReadsForAlignment) {
        delete read;
    }
    tumorReadsForAlignment.clear();
	tumorReadsForAlignmentSize = 0;

    for(pileRead* read : normalReadsForAlignment) {
        delete read;
    }
    normalReadsForAlignment.clear();
	normalReadsForAlignmentSize = 0;

    for(pileRead* read : normalCache) {
        delete read;
    }
    normalCache.clear();

    for(pileRead* read : tumorCache) {
        delete read;
    }
    tumorCache.clear();

    while(!tumorReads.empty())
    {
        delete tumorReads.front();
        tumorReads.pop();
    }

    while(!normalReads.empty())
    {
        delete normalReads.front();
        normalReads.pop();
    }
}