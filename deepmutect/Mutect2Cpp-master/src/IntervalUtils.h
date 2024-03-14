//
// Created by 梦想家xixi on 2021/10/12.
//

#ifndef MUTECT2CPP_MASTER_INTERVALUTILS_H
#define MUTECT2CPP_MASTER_INTERVALUTILS_H

#include "SimpleInterval.h"
#include "samtools/SAMFileHeader.h"

class SimpleInterval;
class Locatable;

class IntervalUtils {
public:
    /**
     * Create a new interval, bounding start and stop by the start and end of contig
     *
     * This function will return null if start and stop cannot be adjusted in any reasonable way
     * to be on the contig.  For example, if start and stop are both past the end of the contig,
     * there's no way to fix this, and null will be returned.
     *
     * @param contig our contig
     * @param start our start as an arbitrary integer (may be negative, etc)
     * @param stop our stop as an arbitrary integer (may be negative, etc)
     * @param contigLength length of the contig
     * @return a valid interval over contig, or null if a meaningful interval cannot be created
     */
    static std::shared_ptr<SimpleInterval> trimIntervalToContig(int contig, int start, int stop, int contigLength);

    /**
    * Tests whether the first Locatable starts after the end of the second Locatable
    *
    * @param first first Locatable
    * @param second second Locatable
    * @param dictionary sequence dictionary used to determine contig ordering
    * @return true if first starts after the end of second, otherwise false
    */
    static bool isAfter(Locatable & first, Locatable & second, SAMSequenceDictionary& dictionary);

    /**
     * Determines the relative contig ordering of first and second using the provided sequence dictionary
     *
     * @param first first Locatable
     * @param second second Locatable
     * @param dictionary sequence dictionary used to determine contig ordering
     * @return 0 if the two contigs are the same, a negative value if first's contig comes before second's contig,
     *         or a positive value if first's contig comes after second's contig
     */
    static int compareContigs(Locatable & first, Locatable & second, SAMSequenceDictionary& dictionary);
};

#endif //MUTECT2CPP_MASTER_INTERVALUTILS_H

