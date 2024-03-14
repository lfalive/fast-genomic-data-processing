//
// The class used to filter out reads for Mutect2Cpp-master
// Created by lhh on 10/19/21.
//

#ifndef MUTECT2CPP_MASTER_READFILTER_H
#define MUTECT2CPP_MASTER_READFILTER_H

#include "htslib/sam.h"
#include "samtools/SAMRecord.h"


class ReadFilter {
public:
    static bool ReadLengthTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool NotSecondaryAlignmentTest(std::shared_ptr<SAMRecord> & originalRead) ;
    static bool GoodCigarTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool NonZeroReferenceLengthAlignmentTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool PassesVendorQualityCheck(std::shared_ptr<SAMRecord> &originalRead);
    static bool MappedTest(std::shared_ptr<SAMRecord> & originalRead );
    static bool MappingQualityAvailableTest(std::shared_ptr<SAMRecord> &originalRead);
    static bool NotDuplicateTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool MappingQualityTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool MappingQualityNotZeroTest(std::shared_ptr<SAMRecord> & originalRead);
    static bool WellformedTest(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header);
    static bool test(std::shared_ptr<SAMRecord> & originalRead, SAMFileHeader* header);

    // overrite some method without using SAMRecord class    // TODO: add NonChimericOriginalAlignmentReadFilter
    static bool ReadLengthTest(bam1_t * read);
    static bool MappingQualityAvailableTest(bam1_t * read);
    static bool MappingQualityTest(bam1_t * read);
    static bool MappingQualityNotZeroTest(bam1_t * read);
    static bool MappedReadTest(bam1_t *read, sam_hdr_t * hdr);
    static bool NotSecondaryAlignmentTest(bam1_t * read);
    static bool NotDuplicateTest(bam1_t * read);
    static bool PassesVendorQualityCheck(bam1_t * read);
    static bool NonZeroReferenceLengthAlignmentTest(bam1_t * read);
    static bool GoodCigarTest(bam1_t * read);
    static bool WellformedTest(bam1_t * read, sam_hdr_t * hdr);
    static bool test(bam1_t * read, sam_hdr_t * hdr);

private:
    static const int MIN_READ_LENGTH = 30;
    static const int READ_QUALITY_FILTER_THRESHOLD = 20;

};


#endif //MUTECT2CPP_MASTER_READFILTER_H
