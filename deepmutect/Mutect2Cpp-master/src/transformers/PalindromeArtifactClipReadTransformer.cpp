//
// Created by 梦想家xixi on 2022/2/22.
//

#include "PalindromeArtifactClipReadTransformer.h"

#include <utility>
#include "read/ReadUtils.h"
#include "utils/BaseUtils.h"
#include "clipping/ReadClipper.h"

const double PalindromeArtifactClipReadTransformer::MIN_FRACTION_OF_MATCHING_BASES = 0.9;

PalindromeArtifactClipReadTransformer::PalindromeArtifactClipReadTransformer(
        std::shared_ptr<ReferenceCache> referenceDataSource, SAMFileHeader *header, int minPalindromeSize) : referenceDataSource(std::move(referenceDataSource)), header(header), minPalindromeSize(minPalindromeSize){

}

std::shared_ptr<SAMRecord> PalindromeArtifactClipReadTransformer::apply(const std::shared_ptr<SAMRecord> &read) {
    int adaptorBoundary = read->getAdaptorBoundary();
    if(!read->isProperlyPaired() || adaptorBoundary == ReadUtils::CANNOT_COMPUTE_ADAPTOR_BOUNDARY) {
        return read;
    }
    std::shared_ptr<Cigar> cigar = read->getCigar();
    CigarOperator firstOperator = cigar->getFirstCigarElement().getOperator();
    CigarOperator lastOperator = cigar->getLastCigarElement().getOperator();
    bool readIsUpstreamOfMate = read->getFragmentLength() > 0;

    if((readIsUpstreamOfMate && firstOperator != S && firstOperator != I) ||
            (!readIsUpstreamOfMate && lastOperator != S && lastOperator != I)) {
        return read;
    }

    int potentialArtifactBaseCount = readIsUpstreamOfMate ? cigar->getFirstCigarElement().getLength() :
                                     cigar->getLastCigarElement().getLength();

    int numBasesToCompare = std::min(potentialArtifactBaseCount + minPalindromeSize, read->getLength());

    int contig = header->getSequenceIndex(read->getContig());
    int refStart = readIsUpstreamOfMate ? adaptorBoundary - numBasesToCompare : adaptorBoundary + 1;
    int refEnd = readIsUpstreamOfMate ? adaptorBoundary - 1 : adaptorBoundary + numBasesToCompare;

    if(refStart < 1 || refEnd > header->getSequenceDictionary().getSequences()[contig].getSequenceLength()) {
        return read;
    }

    if((readIsUpstreamOfMate && refStart < read->getStart()) || (!readIsUpstreamOfMate && read->getEnd() < refEnd)) {
        return read;
    }


    int numMatch = 0;

    int readIndex = readIsUpstreamOfMate ? numBasesToCompare - 1 : read->getLength() - 1;
    int length = 0;
    std::shared_ptr<uint8_t[]> refBases = referenceDataSource->getSubsequenceAt(contig, refStart, refEnd, length);

    for(int i = 0; i < length; i++) {
        if(BaseUtils::getComplement(refBases[i]) == read->getBase(readIndex)) {
            numMatch++;
        }
        readIndex--;
    }

    if(numMatch / static_cast<double>(numBasesToCompare) >= MIN_FRACTION_OF_MATCHING_BASES) {
        ReadClipper readClipper(read);
        ClippingOp clippingOp = readIsUpstreamOfMate ? ClippingOp(0, potentialArtifactBaseCount - 1) :
                ClippingOp(read->getLength() - potentialArtifactBaseCount, read->getLength());
        readClipper.addOp(clippingOp);
        return readClipper.clipRead(HARDCLIP_BASES);
    } else {
        return read;
    }
}

bam1_t * PalindromeArtifactClipReadTransformer::apply(bam1_t * read, sam_hdr_t * hdr)
{
    int adaptorBoundary = ReadUtils::getAdaptorBoundary(read, hdr);
    if(!ReadUtils::isProperlyPaired(read) || adaptorBoundary == ReadUtils::CANNOT_COMPUTE_ADAPTOR_BOUNDARY)
        return read;

    uint32_t * cigar = bam_get_cigar(read);
    uint32_t firstOperator = bam_cigar_op(cigar[0]);
    uint32_t lastOperator = bam_cigar_op(cigar[read->core.n_cigar-1]);
    bool readIsUpstreamOfMate = read->core.isize > 0;

    if((readIsUpstreamOfMate && firstOperator != BAM_CSOFT_CLIP && firstOperator != BAM_CINS) ||
        (!readIsUpstreamOfMate && lastOperator != BAM_CSOFT_CLIP && lastOperator != BAM_CINS)) {
        return read;
    }

    int potentialArtifactBaseCount = readIsUpstreamOfMate ? bam_cigar_oplen(cigar[0]) : bam_cigar_oplen(cigar[read->core.n_cigar-1]);

    int numBasesToCompare = std::min(potentialArtifactBaseCount + minPalindromeSize, read->core.l_qseq);

    int tid = read->core.tid;
    int refStart = readIsUpstreamOfMate ? adaptorBoundary - numBasesToCompare : adaptorBoundary + 1;
    int refEnd = readIsUpstreamOfMate ? adaptorBoundary - 1 : adaptorBoundary + numBasesToCompare;

    if(refStart < 1 || refEnd > sam_hdr_tid2len(hdr, tid)) {
        return read;
    }

    if((readIsUpstreamOfMate && refStart < read->core.pos) || (!readIsUpstreamOfMate && ReadUtils::getEnd(read) < refEnd)) {
        return read;
    }

    int numMatch = 0;

    int readIndex = readIsUpstreamOfMate ? numBasesToCompare - 1 : read->core.l_qseq - 1;
    int length = 0;
    std::shared_ptr<uint8_t[]> refBases = referenceDataSource->getSubsequenceAt(tid, refStart, refEnd, length);

	/*std::string refBasesString((char *) refBases.get(), length);
	std::cout<<refBasesString<<'\n';*/

    uint8_t * bases = bam_get_seq(read);
    for(int i = 0; i < length; i++) {
        if(BaseUtils::getComplement(refBases[i]) == ReadUtils::decodeBase(bam_seqi(bases, readIndex))) {
            numMatch++;
        }
        readIndex--;
    }

    if(numMatch / static_cast<double>(numBasesToCompare) >= MIN_FRACTION_OF_MATCHING_BASES) {

        int start = readIsUpstreamOfMate ? 0 : read->core.l_qseq - potentialArtifactBaseCount;
        int stop = readIsUpstreamOfMate ? potentialArtifactBaseCount - 1 : read->core.l_qseq;
        return ReadClipper::hardClipRead(read, start, stop, hdr);

    } else {
        return read;
    }

}
