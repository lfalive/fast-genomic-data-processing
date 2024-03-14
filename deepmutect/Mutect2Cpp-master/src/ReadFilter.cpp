//
// Created by lhh on 10/19/21.
//

#include "ReadFilter.h"
#include "read/CigarUtils.h"
#include "QualityUtils.h"
#include "read/ReadUtils.h"


bool ReadFilter::NotSecondaryAlignmentTest(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->isSecondaryAlignment();
}

bool ReadFilter::GoodCigarTest(std::shared_ptr<SAMRecord> & originalRead) {
    return CigarUtils::isGood(originalRead->getCigar());
}

bool ReadFilter::GoodCigarTest(bam1_t * read)
{
    return CigarUtils::isGood(read);
}

bool ReadFilter::NonZeroReferenceLengthAlignmentTest(std::shared_ptr<SAMRecord> & originalRead) {
    for(const CigarElement& element : originalRead->getCigarElements()) {
        if(CigarOperatorUtils::getConsumesReferenceBases(element.getOperator()) && element.getLength() > 0) {
            return true;
        }
    }
    return false;
}

bool ReadFilter::NonZeroReferenceLengthAlignmentTest(bam1_t * read) {
    uint32_t * cigarArray = bam_get_cigar(read);
    for(int i=0; i<read->core.n_cigar; i++)
    {
        uint32_t cigarElement = cigarArray[i];
        if(ReadUtils::consumesReferenceBases(bam_cigar_op(cigarElement)) && bam_cigar_oplen(cigarElement) > 0)
            return true;
    }
    return false;
}

bool ReadFilter::PassesVendorQualityCheck(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->failsVendorQualityCheck();
}

bool ReadFilter::MappedTest(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->isUnmapped();
}

bool ReadFilter::MappingQualityAvailableTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getMappingQuality() != QualityUtils::MAPPING_QUALITY_UNAVALIABLE;
}

bool ReadFilter::MappingQualityAvailableTest(bam1_t * read){
    return read->core.qual != QualityUtils::MAPPING_QUALITY_UNAVALIABLE;
}

bool ReadFilter::NotDuplicateTest(std::shared_ptr<SAMRecord> & originalRead) {
    return ! originalRead->isDuplicate();
}

bool ReadFilter::MappingQualityTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getMappingQuality() >= READ_QUALITY_FILTER_THRESHOLD;
}

bool ReadFilter::MappingQualityTest(bam1_t * read) {
    return read->core.qual >= READ_QUALITY_FILTER_THRESHOLD;
}

bool ReadFilter::MappingQualityNotZeroTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getMappingQuality() != 0;
}

bool ReadFilter::WellformedTest(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header) {
    return (originalRead->isUnmapped() || originalRead->getStart() > 0) &&
            (originalRead->isUnmapped() || (originalRead->getEnd() - originalRead->getStart() + 1) >= 0) &&
            ReadUtils::alignmentAgreesWithHeader(header, originalRead) &&
            // ! originalRead.getReadGroup().empty() &&
            originalRead->getLength() == originalRead->getBaseQualitiesLength() &&
            (originalRead->isUnmapped() || originalRead->getLength() == Cigar::getReadLength(originalRead->getCigarElements())) &&
            (originalRead->getLength() > 0) &&
            (! CigarUtils::containsNOperator(originalRead->getCigarElements()));
}

bool ReadFilter::WellformedTest(bam1_t * read, sam_hdr_t * hdr)
{
    return (ReadUtils::isUnmapped(read, hdr) || read->core.pos + 1 > 0) &&
    (ReadUtils::isUnmapped(read, hdr) || (ReadUtils::getEnd(read) - read->core.pos + 1 >= 0) ) &&
    ReadUtils::alignmentAgreesWithHeader(hdr, read) &&
    read->core.l_qseq == bam_cigar2qlen(read->core.n_cigar, bam_get_cigar(read)) &&
    (ReadUtils::isUnmapped(read, hdr) || read->core.l_qseq == Cigar::getReadLength(read->core.n_cigar, bam_get_cigar(read))) &&
    (read->core.l_qseq > 0) &&
    (!CigarUtils::containsNOperator(read->core.n_cigar, bam_get_cigar(read)));
}

bool ReadFilter::test(bam1_t * read, sam_hdr_t * hdr)
{
//	bool ret1 = ReadLengthTest(read);
//	bool ret2 = NonZeroReferenceLengthAlignmentTest(read);
//	bool ret3 = NotDuplicateTest(read);
//	bool ret4 = NotSecondaryAlignmentTest(read);
//	bool ret5 = GoodCigarTest(read);
//	bool ret6 = PassesVendorQualityCheck(read);
//	bool ret7 = MappedReadTest(read, hdr);
//	bool ret8 = MappingQualityAvailableTest(read);
//	bool ret9 = MappingQualityNotZeroTest(read);
//	bool ret10 = MappingQualityTest(read);
//	bool ret11 = WellformedTest(read, hdr);
    return ReadLengthTest(read) && NonZeroReferenceLengthAlignmentTest(read) && NotDuplicateTest(read) &&
            NotSecondaryAlignmentTest(read) && GoodCigarTest(read) && PassesVendorQualityCheck(read) && MappedReadTest(read, hdr) &&
            MappingQualityAvailableTest(read) && MappingQualityNotZeroTest(read) && MappingQualityTest(read) &&
            WellformedTest(read, hdr);
}

bool ReadFilter::test(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header) {
//    bool ret1 = ReadLengthTest();
//    bool ret2 = NonZeroReferenceLengthAlignmentTest();
//    bool ret3 = NotDuplicateTest();
//    bool ret4 = NotSecondaryAlignmentTest();
//    bool ret5 = GoodCigarTest();
//    bool ret6 = PassesVendorQualityCheck();
//    bool ret7 = MappedTest();
//    bool ret8 = MappingQualityAvailableTest();
//    bool ret9 = MappingQualityNotZeroTest();
//    bool ret10 = MappingQualityTest();
//    bool ret11 = WellformedTest();
    return ReadLengthTest(originalRead)&& NonZeroReferenceLengthAlignmentTest(originalRead) && NotDuplicateTest(originalRead) && NotSecondaryAlignmentTest(originalRead) &&
    GoodCigarTest(originalRead) && PassesVendorQualityCheck(originalRead) && MappedTest(originalRead) && MappingQualityAvailableTest(originalRead) && MappingQualityNotZeroTest(originalRead) && MappingQualityTest(originalRead) && WellformedTest(originalRead,header);

//    return true;
}

bool ReadFilter::ReadLengthTest(std::shared_ptr<SAMRecord> & originalRead) {
    return originalRead->getLength() >= MIN_READ_LENGTH && originalRead->getLength() < INT32_MAX;
}

bool ReadFilter::ReadLengthTest(bam1_t *read) {
    return read->core.l_qseq >= MIN_READ_LENGTH && read->core.l_qseq < INT32_MAX;
}

bool ReadFilter::MappingQualityNotZeroTest(bam1_t *read)
{
    return (read->core.qual != 0);
}


bool ReadFilter::MappedReadTest(bam1_t *read, sam_hdr_t * hdr)
{
    return !ReadUtils::isUnmapped(read, hdr);
}

bool ReadFilter::NotSecondaryAlignmentTest(bam1_t *read)
{
    return (read->core.flag & BAM_FSECONDARY) == 0;
}

bool ReadFilter::NotDuplicateTest(bam1_t *read)
{
    return (read->core.flag & BAM_FDUP) == 0;
}

bool ReadFilter::PassesVendorQualityCheck(bam1_t *read)
{
    return (read->core.flag & BAM_FQCFAIL) == 0;
}

