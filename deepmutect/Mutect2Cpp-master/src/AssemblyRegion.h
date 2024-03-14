//
// Created by 梦想家xixi on 2021/10/14.
//

#ifndef MUTECT2CPP_MASTER_ASSEMBLYREGION_H
#define MUTECT2CPP_MASTER_ASSEMBLYREGION_H

#include <vector>
#include <set>
#include <htslib/sam.h>
#include "SimpleInterval.h"
#include "ActivityProfileState.h"
#include "samtools/SAMRecord.h"
#include "ReferenceCache.h"

class AssemblyRegion : public Locatable{
private:
    SAMFileHeader *header = nullptr;

    /**
     * The reads included in this assembly region.  May be empty upon creation, and expand / contract
     * as reads are added or removed from this region.
     */
    std::vector<std::shared_ptr<SAMRecord>> reads;

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this assembly region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     */
     std::vector<std::shared_ptr<ActivityProfileState>> supportingStates;

    /**
    * The raw span of this assembly region, not including the region extension
    */
    std::shared_ptr<SimpleInterval> activeRegionLoc;

    /**
     * The span of this assembly region on the genome, including the region extension
     */
     std::shared_ptr<SimpleInterval> extendedLoc;

    /**
    * The extension, in bp, of this region. The extension is >= 0 bp in size, and indicates how much padding was
    * requested for the region.
    */
    const int extension;

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
     bool isActive;

    /**
    * The span of this assembly region, including the bp covered by all reads in this
    * region.  This union of extensionLoc and the loc of all reads in this region.
    *
    * Must be at least as large as extendedLoc, but may be larger when reads
    * partially overlap this region.
    */
    std::shared_ptr<SimpleInterval> spanIncludingReads;

    /**
     * Indicates whether the region has been finalized
     */
     bool hasBeenFinalized = false;

    /**
    * Create a new SimpleInterval, bounding start and stop by the start and end of contig
    *
    * This function will return null if start and stop cannot be adjusted in any reasonable way
    * to be on the contig.  For example, if start and stop are both past the end of the contig,
    * there's no way to fix this, and null will be returned.
    *
    * @param contig our contig
    * @param start our start as an arbitrary integer (may be negative, etc)
    * @param stop our stop as an arbitrary integer (may be negative, etc)
    * @return a valid genome loc over contig, or null if a meaningful genome loc cannot be created
    */
    std::shared_ptr<SimpleInterval> trimIntervalToContig(int contig, int start, int stop);

    void checkStates(SimpleInterval& activeRegion);

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases.
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @param genomeLoc a non-null genome loc indicating the base span of the bp we'd like to get the reference for
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    std::shared_ptr<uint8_t[]> getReference(ReferenceCache* referenceReader, int padding, const std::shared_ptr<SimpleInterval> & genomeLoc, int & length);

public:
    AssemblyRegion(SimpleInterval const &activeRegionLoc, std::vector<std::shared_ptr<ActivityProfileState>>  supportingStates, bool isActive, int extension, SAMFileHeader * header);

    /**
     * Simple interface to create an assembly region that isActive without any profile state
     */
    AssemblyRegion(SimpleInterval const &activeRegionLoc, int extension);
    ~AssemblyRegion();

    std::string getContig() const override {return activeRegionLoc->getContig();}

    int getStart() const override {return activeRegionLoc->getStart();}

    int getEnd() const override {return activeRegionLoc->getEnd();}

    friend std::ostream & operator<<(std::ostream & os, AssemblyRegion & assemblyRegion);

    /**
     * Does this region represent an active region (all isActiveProbs above threshold) or
     * an inactive region (all isActiveProbs below threshold)?
     */
    bool getIsActive() const {return isActive;}

    /**
     * Override activity state of the region
     *
     * Note: Changing the isActive state after construction is a debug-level operation that only engine classes
     * like AssemblyRegionWalker should be able to do
     *
     * @param value new activity state of this region
     */
    void setIsActive(bool value);

    /**
     * Get the span of this assembly region including the extension value
     * @return a non-null SimpleInterval
     */
    std::shared_ptr<SimpleInterval> getExtendedSpan() {return extendedLoc;}

    /**
     * Get the raw span of this assembly region (excluding the extension)
     * @return a non-null SimpleInterval
     */
    std::shared_ptr<SimpleInterval> getSpan() {return activeRegionLoc;}

    /**
    * Get an unmodifiable copy of the list of reads currently in this assembly region.
    *
    * The reads are sorted by their coordinate position.
    * @return an unmodifiable and inmutable copy of the reads in the assembly region.
   */
    std::vector<std::shared_ptr<SAMRecord>> & getReads();

    /**
     * Returns the header for the reads in this region.
     */
    SAMFileHeader* getHeader();

    /**
     * Intersect this assembly region with the allowed intervals, returning a list of active regions
     * that only contain locations present in intervals
     *
     * Note: modifications to the returned list have no effect on this region object.
     *
     * Note that the returned list may be empty, if this active region doesn't overlap the set at all
     *
     * Note that the resulting regions are all empty, regardless of whether the current active region has reads
     *
     * @param intervals a non-null set of intervals that are allowed
     * @return an ordered list of active region where each interval is contained within intervals
     */
     //TODO:std::vector<AssemblyRegion> splitAndTrimToIntervals(std::set<SimpleInterval>);

    /**
    * Trim this region to just the span, producing a new assembly region without any reads that has only
    * the extent of newExtend intersected with the current extent
    * @param span the new extend of the active region we want
    * @param extensionSize the extensionSize size we want for the newly trimmed active region
    * @return a non-null, empty assembly region
    */
    std::shared_ptr<AssemblyRegion> trim(const std::shared_ptr<SimpleInterval>& span, const std::shared_ptr<SimpleInterval>& extendedSpan);

    /**
     * Get the extension applied to this region
     *
     * The extension is >= 0 bp in size, and indicates how much padding was requested for the region
     *
     * @return the size in bp of the region extension
     */
     int getExtension() const { return extension; }

     /*TODO:bool readOverlapsRegion( GATKRead read);
        void addAll
        byte[] getFullReference
        static byte[] getReference
        byte[] getAssemblyRegionReference
        byte[] getAssemblyRegionReference */


    /**
    * The span of this assembly region, including the bp covered by all reads in this
    * region.  This union of extensionLoc and the loc of all reads in this region.
    *
    * Must be at least as large as extendedLoc, but may be larger when reads
    * partially overlap this region.
    */
    std::shared_ptr<SimpleInterval> getReadSpanLoc() const {return spanIncludingReads;}

    /**
     * An ordered list (by genomic coordinate) of the ActivityProfileStates that went
     * into this active region.  May be empty, which says that no supporting states were
     * provided when this region was created.
     * The returned list is unmodifiable.
     */
    std::vector<std::shared_ptr<ActivityProfileState>> getSupportingStates() {return supportingStates;}

    /**
    * Is this region equal to other, excluding any reads in either region in the comparison
    * @param other the other active region we want to test
    * @return true if this region is equal, excluding any reads and derived values, to other
    */
    bool equalsIgnoreReads(AssemblyRegion const &other);

    void setFinalized(bool value);

    bool isFinalized() const {return hasBeenFinalized;}

    void setRead(std::vector<std::shared_ptr<SAMRecord>> & reads);

    bool readOverlapsRegion(std::shared_ptr<SAMRecord> & read);

    /**
     * Clear all of the reads currently in this region
     */
    void clearReads();

    /**
     * Add all readsToAdd to this region
     * @param readsToAdd a collection of readsToAdd to add to this active region
     */
    void addAll(std::vector<std::shared_ptr<SAMRecord>> & readsToAdd);

    /**
     * Get the reference bases from referenceReader spanned by the extended location of this active region,
     * including additional padding bp on either side.  If this expanded region would exceed the boundaries
     * of the active region's contig, the returned result will be truncated to only include on-genome reference
     * bases
     *
     * @param referenceReader the source of the reference genome bases
     * @param padding the padding, in BP, we want to add to either side of this active region extended region
     * @return a non-null array of bytes holding the reference bases in referenceReader
     */
    std::shared_ptr<uint8_t[]> getAssemblyRegionReference(ReferenceCache * cache, int padding, int & length);

    /**
     * Remove all of the reads in readsToRemove from this region
     * @param readsToRemove the set of reads we want to remove
     */
    void removeAll(const std::vector<std::shared_ptr<SAMRecord>>& readsToRemove);

    /**
    * Add read to this region
    *
    * Read must have alignment start >= than the last read currently in this active region.
    *
    * @throws IllegalArgumentException if read doesn't overlap the extended region of this active region
    *
    * @param read a non-null GATKRead
    */
    void add(std::shared_ptr<SAMRecord>& read);

	void printRegionInfo();

	void sortReadsByCoordinate();

    std::shared_ptr<uint8_t[]>  getFullReference(ReferenceCache * cache, int padding, int & length);
};


#endif //MUTECT2CPP_MASTER_ASSEMBLYREGION_H
