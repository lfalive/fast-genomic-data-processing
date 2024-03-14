/**
 * The implementation of ReferenceCache class
 */
#include <iostream>
#include "ReferenceCache.h"


ReferenceCache::ReferenceCache(char * refName, SAMFileHeader* header, int tid) : tid(tid), header(header)
{
    fai = fai_load3_format(refName, nullptr, nullptr, FAI_CREATE, FAI_FASTA);
    start = 0;
    end = header->getSequenceDictionary().getSequences()[tid].getSequenceLength();  // TODO: make it more elegan
    std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start + 1) + '-' + std::to_string(end + 1);
    bases = fai_fetch64(fai, region.c_str(), &len);
	for (int i = 0; i < end; ++i) {
		bases[i] = (char)std::toupper(bases[i]);   // toUpperCase
		if (bases[i]!='A' && bases[i]!='C' && bases[i]!='G' && bases[i]!='T')   // convertIUPACtoN
			bases[i] = 'N';
	}
}

ReferenceCache::~ReferenceCache()
{
    clear();
    fai_destroy(fai);
}

std::string ReferenceCache::getContig()
{
    return header->getSequenceDictionary().getSequences()[tid].getSequenceName();
}

void ReferenceCache::clear()
{
    if (this->bases && strlen(this->bases)) // strlen does not check for null pointer
        free(bases);
}




char ReferenceCache::getBase(hts_pos_t pos)
{
/*    while(pos > end) {
        advanceLoad();
    }*/
    /*char base = bases[pos - start];
    if(base >= 'A' && base <= 'Z')
        return base;
    else if(base >= 'a' && base <= 'z')
        return base - ('a' - 'A');*/
	return bases[pos - start];
}

void ReferenceCache::advanceLoad() {
    start = end + 1;
    end = std::min(start + 999999, static_cast<hts_pos_t>(header->getSequenceDictionary().getSequences()[tid].getSequenceLength()-1));
    std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start + 1) + '-' + std::to_string(end + 1);
    hts_pos_t seq_len;
    clear();
    bases = fai_fetch64(fai, region.c_str(), &len);
}

void ReferenceCache::setTid(int tid) {
    this -> tid = tid;
    start = 0;
    end = std::min(99999, header->getSequenceDictionary().getSequences()[tid].getSequenceLength()-1);
    std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start+1) + '-' + std::to_string(end+1);
    clear();
    bases = fai_fetch64(fai, region.c_str(), &len);

}

std::shared_ptr<uint8_t[]> ReferenceCache::getSubsequenceAt(int tid, int start, int stop, int & length) {
    if(tid == this->tid && start >= this->start && stop <= this->end) {
        std::shared_ptr<uint8_t[]> ret(new uint8_t[stop-start+1]{0});
        std::copy(bases+start-this->start, bases+stop-this->start + 1, ret.get());
        length = stop - start + 1;
        return ret;
    }
    else {
        std::string region = header->getSequenceDictionary().getSequences()[tid].getSequenceName() + ':' + std::to_string(start+1) + '-' + std::to_string(stop+1);
        hts_pos_t seq_len;
        auto * ret = reinterpret_cast<uint8_t*>(fai_fetch64(fai, region.c_str(), &seq_len));
        std::shared_ptr<uint8_t[]> toRet(new uint8_t[seq_len]);
        std::copy(ret, ret + seq_len, toRet.get());
        length = seq_len;
        free(ret);
        return toRet;
    }
}
