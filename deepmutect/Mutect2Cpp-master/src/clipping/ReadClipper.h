//
// Created by 梦想家xixi on 2021/12/18.
//

#ifndef MUTECT2CPP_MASTER_READCLIPPER_H
#define MUTECT2CPP_MASTER_READCLIPPER_H

#include "samtools/SAMRecord.h"
#include "ClippingOp.h"
#include "ClippingRepresentation.h"

class ClippingOp;

class cigarShift {
public:
	std::vector<uint32_t> *cigar;
	int start;
	int end;

	cigarShift(std::vector<uint32_t> *cigar, int start, int end) {
		this->cigar = cigar;
		this->start = start;
		this->end = end;
	}
};


class ReadClipper {
public:
	std::shared_ptr<SAMRecord> read;
	bool wasClipped;
	std::vector<ClippingOp> ops;

	explicit ReadClipper(std::shared_ptr<SAMRecord> read);

	void addOp(const ClippingOp &op);

	static std::shared_ptr<SAMRecord>
	hardClipToRegion(const std::shared_ptr<SAMRecord> &read, int refStart, int refStop);

	static std::shared_ptr<SAMRecord>
	hardClipBothEndsByReferenceCoordinates(std::shared_ptr<SAMRecord> read, int left, int right);

	static std::shared_ptr<SAMRecord>
	hardClipByReferenceCoordinatesLeftTail(std::shared_ptr<SAMRecord> read, int refStop);

	static std::shared_ptr<SAMRecord>
	hardClipByReferenceCoordinatesRightTail(std::shared_ptr<SAMRecord> read, int refStop);

	static std::shared_ptr<SAMRecord> hardClipLowQualEnds(std::shared_ptr<SAMRecord> read, uint8_t lowQual);

	static std::shared_ptr<SAMRecord> hardClipSoftClippedBases(std::shared_ptr<SAMRecord> read);

	static std::shared_ptr<SAMRecord> revertSoftClippedBases(std::shared_ptr<SAMRecord> read);

	static std::shared_ptr<SAMRecord> hardClipAdaptorSequence(std::shared_ptr<SAMRecord> read);

	std::shared_ptr<SAMRecord> clipRead(ClippingRepresentation algorithm);

	// clip reads, overrite some method in ReadClipper class
	static bam1_t *hardClipRead(bam1_t *read, int start, int stop, sam_hdr_t *hdr);

	static bam1_t *applyHardClipBases(bam1_t *read, int start, int stop, sam_hdr_t *hdr);

	static cigarShift hardClipCigar(uint32_t *cigar, int n_cigar, int start, int stop);

	/**
	  * Checks if a hard clipped cigar left a read starting or ending with deletions or gap (N)
	  * and cleans it up accordingly.
	  *
	  * @param cigar the original cigar
	  * @return an object with the shifts (see CigarShift class)
	  */
	static cigarShift cleanHardClippedCigar(std::vector<uint32_t> &newCigar);

private:
	static std::shared_ptr<SAMRecord>
	hardClipToRegion(std::shared_ptr<SAMRecord> read, int refStart, int refStop, int alignmentStart, int alignmentStop);

	std::shared_ptr<SAMRecord> hardClipBothEndsByReferenceCoordinates(int left, int right);

	std::shared_ptr<SAMRecord>
	clipByReferenceCoordinates(int refStart, int refStop, ClippingRepresentation clippingOp, bool runAsserts);

	std::shared_ptr<SAMRecord> clipRead(ClippingRepresentation algorithm, bool runAsserts);

	std::shared_ptr<SAMRecord> hardClipByReferenceCoordinatesLeftTail(int refStop);

	std::shared_ptr<SAMRecord> hardClipLowQualEnds(uint8_t lowQual);

	std::shared_ptr<SAMRecord> clipLowQualEnds(ClippingRepresentation algorithm, uint8_t lowQual);

	std::shared_ptr<SAMRecord> hardClipSoftClippedBases();

	std::shared_ptr<SAMRecord> revertSoftClippedBases();

	std::shared_ptr<SAMRecord> hardClipAdaptorSequence();


};


#endif //MUTECT2CPP_MASTER_READCLIPPER_H
