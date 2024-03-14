//
// Created by 梦想家xixi on 2021/12/18.
//

#include "ReadClipper.h"
#include <utility>
#include "read/ReadUtils.h"

std::shared_ptr<SAMRecord>
ReadClipper::hardClipToRegion(const std::shared_ptr<SAMRecord> &read, int refStart, int refStop) {
	int start = read->getStart();
	int stop = read->getEnd();
	return hardClipToRegion(read, refStart, refStop, start, stop);
}

std::shared_ptr<SAMRecord>
ReadClipper::hardClipToRegion(std::shared_ptr<SAMRecord> read, int refStart, int refStop, int alignmentStart,
                              int alignmentStop) {
	if (alignmentStart <= refStop && alignmentStop >= refStart) {
		if (alignmentStart < refStart && alignmentStop > refStop)
			return hardClipBothEndsByReferenceCoordinates(read, refStart - 1, refStop + 1);

		if (alignmentStart < refStart)
			return hardClipByReferenceCoordinatesLeftTail(read, refStart - 1);

		if (alignmentStop > refStop)
			return hardClipByReferenceCoordinatesRightTail(read, refStop + 1);

		return read;
	}
	return ReadUtils::emptyRead(read);
}

ReadClipper::ReadClipper(std::shared_ptr<SAMRecord> read) : read(std::move(read)), wasClipped(false) {
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipBothEndsByReferenceCoordinates(const int left, const int right) {
	if (read->getLength() == 0 || left == right) {
		return ReadUtils::emptyRead(read);
	}
	std::shared_ptr<SAMRecord> leftTailRead = clipByReferenceCoordinates(right, -1, HARDCLIP_BASES, true);
	if (left > leftTailRead->getEnd()) {
		return ReadUtils::emptyRead(read);
	}
	ReadClipper clipper = ReadClipper(leftTailRead);
	return clipper.hardClipByReferenceCoordinatesLeftTail(left);
}

std::shared_ptr<SAMRecord>
ReadClipper::clipByReferenceCoordinates(int refStart, int refStop, ClippingRepresentation clippingOp, bool runAsserts) {
	if (read->getLength() == 0) {
		return read;
	}
	if (clippingOp == SOFTCLIP_BASES && read->isUnmapped()) {
		throw std::invalid_argument("Cannot soft-clip read by reference coordinates because it is unmapped");
	}
	int start, stop;

	if (refStart < 0) {
		if (refStop < 0)
			throw std::invalid_argument("Only one of refStart or refStop must be < 0");

		start = 0;
		stop = ReadUtils::getReadCoordinateForReferenceCoordinate(read, refStop, LEFT_TAIL);
	} else {
		if (refStop >= 0)
			throw std::invalid_argument("Either refStart or refStop must be < 0");

		start = ReadUtils::getReadCoordinateForReferenceCoordinate(read, refStart, RIGHT_TAIL);
		stop = read->getLength() - 1;
	}

	if (start < 0 || stop > read->getLength() - 1)
		throw std::invalid_argument("Trying to clip before the start or after the end of a read");

	if (start > stop)
		throw std::invalid_argument("START > STOP -- this should never happen, please check read");

	if (start > 0 && stop < read->getLength() - 1)
		throw std::invalid_argument("Trying to clip the middle of the read");

	addOp(ClippingOp(start, stop));
	std::shared_ptr<SAMRecord> clippedRead = clipRead(clippingOp, runAsserts);
	ops.clear();
	return clippedRead;
}

void ReadClipper::addOp(const ClippingOp &op) {
	ops.emplace_back(op);
}

std::shared_ptr<SAMRecord> ReadClipper::clipRead(ClippingRepresentation algorithm, bool runAsserts) {
	Mutect2Utils::validateArg(algorithm != NULL_ClippingRepresentation, "null is not allowed there.");
	if (ops.empty())
		return read;

	std::shared_ptr<SAMRecord> clippedRead = read;
	for (ClippingOp op: ops) {
		int readLength = clippedRead->getLength();
		if (op.start < readLength) {
			ClippingOp fixedOperation = op;
			if (op.stop >= readLength) {
				fixedOperation = ClippingOp(op.start, readLength - 1);
			}
			clippedRead = fixedOperation.apply(algorithm, clippedRead, runAsserts);
		}
	}
	wasClipped = true;
	ops.clear();
	if (clippedRead->isEmpty())
		return ReadUtils::emptyRead(clippedRead);
	return clippedRead;
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipByReferenceCoordinatesLeftTail(int refStop) {
	return clipByReferenceCoordinates(-1, refStop, HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord>
ReadClipper::hardClipBothEndsByReferenceCoordinates(std::shared_ptr<SAMRecord> read, int left, int right) {
	return ReadClipper(std::move(read)).hardClipBothEndsByReferenceCoordinates(left, right);
}

std::shared_ptr<SAMRecord>
ReadClipper::hardClipByReferenceCoordinatesLeftTail(std::shared_ptr<SAMRecord> read, int refStop) {
	return ReadClipper(std::move(read)).clipByReferenceCoordinates(-1, refStop, HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord>
ReadClipper::hardClipByReferenceCoordinatesRightTail(std::shared_ptr<SAMRecord> read, int refStart) {
	return ReadClipper(std::move(read)).clipByReferenceCoordinates(refStart, -1, HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::clipLowQualEnds(ClippingRepresentation algorithm, uint8_t lowQual) {
	if (read->isEmpty())
		return read;

	int readLength = read->getLength();
	int leftClipIndex = 0;
	int rightClipIndex = readLength - 1;

	while (rightClipIndex >= 0 && read->getBaseQuality(rightClipIndex) <= lowQual) {
		rightClipIndex--;
	}
	while (leftClipIndex < readLength && read->getBaseQuality(leftClipIndex) <= lowQual) {
		leftClipIndex++;
	}
	if (leftClipIndex > rightClipIndex)
		return ReadUtils::emptyRead(read);

	if (rightClipIndex < readLength - 1) {
		addOp(ClippingOp(rightClipIndex + 1, readLength - 1));
	}

	if (leftClipIndex > 0) {
		addOp(ClippingOp(0, leftClipIndex - 1));
	}
	return clipRead(algorithm, true);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipLowQualEnds(uint8_t lowQual) {
	return clipLowQualEnds(HARDCLIP_BASES, lowQual);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipLowQualEnds(std::shared_ptr<SAMRecord> read, uint8_t lowQual) {
	return ReadClipper(std::move(read)).hardClipLowQualEnds(lowQual);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipSoftClippedBases(std::shared_ptr<SAMRecord> read) {
	return ReadClipper(std::move(read)).hardClipSoftClippedBases();
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipSoftClippedBases() {
	if (read->isEmpty())
		return read;

	int readIndex = 0;
	int cutLeft = -1;            // first position to hard clip (inclusive)
	int cutRight = -1;           // first position to hard clip (inclusive)
	bool rightTail = false;

	for (CigarElement cigarElement: read->getCigarElements()) {
		if (cigarElement.getOperator() == S) {
			if (rightTail) {
				cutRight = readIndex;
			} else {
				cutLeft = readIndex + cigarElement.getLength() - 1;
			}
		} else if (cigarElement.getOperator() != H) {
			rightTail = true;
		}

		if (CigarOperatorUtils::getConsumesReadBases(cigarElement.getOperator())) {
			readIndex += cigarElement.getLength();
		}
	}

	if (cutRight >= 0) {
		addOp(ClippingOp(cutRight, read->getLength() - 1));
	}
	if (cutLeft >= 0) {
		addOp(ClippingOp(0, cutLeft));
	}
	return clipRead(HARDCLIP_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::revertSoftClippedBases(std::shared_ptr<SAMRecord> read) {
	return ReadClipper(std::move(read)).revertSoftClippedBases();
}

std::shared_ptr<SAMRecord> ReadClipper::revertSoftClippedBases() {
	if (read->isEmpty())
		return read;
	addOp(ClippingOp(0, 0));
	return clipRead(REVERT_SOFTCLIPPED_BASES, true);
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipAdaptorSequence(std::shared_ptr<SAMRecord> read) {
	return ReadClipper(std::move(read)).hardClipAdaptorSequence();
}

std::shared_ptr<SAMRecord> ReadClipper::hardClipAdaptorSequence() {
	int adaptorBoundary = read->getAdaptorBoundary();

	if (adaptorBoundary == ReadUtils::CANNOT_COMPUTE_ADAPTOR_BOUNDARY ||
	    !ReadUtils::isInsideRead(read, adaptorBoundary)) {
		return read;
	}

	return read->isReverseStrand() ? hardClipByReferenceCoordinatesLeftTail(adaptorBoundary)
	                               : hardClipByReferenceCoordinatesRightTail(read, adaptorBoundary);
}

std::shared_ptr<SAMRecord> ReadClipper::clipRead(ClippingRepresentation algorithm) {
	return clipRead(algorithm, true);
}


cigarShift ReadClipper::hardClipCigar(uint32_t *cigar, int n_cigar, int start, int stop) {
	//uint32_t * newCigar = NULL;
	std::vector<uint32_t> newCigar;
	int index = 0;
	int totalHardClipCount = stop - start + 1;
	int i = 0;

	// hard clip the beginning of the cigar string
	if (start == 0) {
		while (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
			totalHardClipCount += bam_cigar_oplen(cigar[i]);
			if (i < n_cigar - 1)
				i++;
			else
				throw std::invalid_argument(
						"Read is entirely hard-clipped, shouldn't be trying to clip it's cigar string");
		}

		// keep clipping until we hit stop
		while (index <= stop) {
			int shift = 0;
			if (ReadUtils::consumesReadBases(cigar[i]))
				shift = bam_cigar_oplen(cigar[i]);

			// we're still clipping or just finished perfectly
			if (index + shift == stop + 1) {
				newCigar.push_back(bam_cigar_gen(totalHardClipCount, BAM_CHARD_CLIP));
			} else if (index + shift > stop + 1) { // element goes beyond what we need to clip
				int elementLengthAfterChopping = bam_cigar_oplen(cigar[i]) - (stop - index + 1);
				newCigar.push_back(bam_cigar_gen(totalHardClipCount, BAM_CHARD_CLIP));
				newCigar.push_back(bam_cigar_gen(elementLengthAfterChopping, bam_cigar_op(cigar[i])));
			}
			index += shift;

			if (index <= stop && i < n_cigar - 1) {
				i++;
			} else {
				break;
			}
		}

		// add the remaining cigar elements
		while (i < n_cigar - 1) {
			i++;
			newCigar.push_back(cigar[i]);
		}
	} else {   // hard clip the end of the cigar string
		// Keep marching on until we find the start
		while (index < start) {
			int shift = 0;
			if (ReadUtils::consumesReadBases(cigar[i]))
				shift = bam_cigar_oplen(cigar[i]);

			// we haven't gotten to the start yet, keep everything as is.
			if (index + shift < start) {
				newCigar.push_back(cigar[i]);
			}// element goes beyond our clip starting position
			else {
				int elementLengthAfterChopping = start - index;
				// if this last element is a HARD CLIP operator, just merge it with our hard clip operator to be added later
				if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
					totalHardClipCount += elementLengthAfterChopping;
				}// otherwise, maintain what's left of this last operator
				else {
					newCigar.push_back(bam_cigar_gen(elementLengthAfterChopping, bam_cigar_op(cigar[i])));
				}
			}

			index += shift;
			if (index < start && i < n_cigar - 1) {
				i++;
			} else {
				break;
			}
		}

		// check if we are hard clipping indels
		while (i < n_cigar - 1) {
			i++;
			if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
				totalHardClipCount += bam_cigar_oplen(cigar[i]);
			}
		}

		newCigar.push_back(bam_cigar_gen(totalHardClipCount, BAM_CHARD_CLIP));

	}
	return cleanHardClippedCigar(newCigar);
}

cigarShift ReadClipper::cleanHardClippedCigar(std::vector<uint32_t> &newCigar) {
	auto *cleanCigar = new std::vector<uint32_t>;
	int shiftFromStart = 0;
	int shiftFromEnd = 0;
	std::stack<uint32_t> cigarStack, inverseCigarStack;

	for (uint32_t cigarElement: newCigar) {
		cigarStack.push(cigarElement);
	}

	for (int passInt = FIRST; passInt != END; passInt++) {
		Passes pass = static_cast<Passes>(passInt);

		int shift = 0;
		int totalHardClip = 0;
		bool readHasStarted = false;
		bool addedHardClips = false;

		while (!cigarStack.empty()) {
			uint32_t cigarElement = cigarStack.top();
			cigarStack.pop();

			if (!readHasStarted && bam_cigar_op(cigarElement) == BAM_CHARD_CLIP) {
				totalHardClip += bam_cigar_oplen(cigarElement);
			}

			// Deletions (D) and gaps (N) are not hardclips (H) and do not consume read bases....
			// so they gets dropped from the edges of the read since readHasStarted is still false.

			readHasStarted |= ReadUtils::consumesReadBases(cigarElement);

			if (readHasStarted) {
				switch (pass) {
					case FIRST:
						if (!addedHardClips && totalHardClip > 0) {
							inverseCigarStack.push(bam_cigar_gen(totalHardClip, BAM_CHARD_CLIP));
						}
						inverseCigarStack.push(cigarElement);
						break;
					case SECOND:
						if (!addedHardClips && totalHardClip > 0) {
							cleanCigar->push_back(bam_cigar_gen(totalHardClip, BAM_CHARD_CLIP));
						}
						cleanCigar->push_back(cigarElement);
						break;
					case END:
						break;
				}
				addedHardClips = true;
			}
		}

		switch (pass) {
			// first pass is from end to start of the cigar elements
			case FIRST:
				shiftFromEnd = shift;
				cigarStack = inverseCigarStack;
				break;
			case SECOND:
				// second pass is from start to end with the end already cleaned
				shiftFromStart = shift;
				break;
			case END:
				break;
		}
	}

	return cigarShift{cleanCigar, shiftFromStart, shiftFromEnd};
}

bam1_t *ReadClipper::applyHardClipBases(bam1_t *read, int start, int stop, sam_hdr_t *hdr) {

	// If the read is unmapped there is no Cigar string and neither should we create a new cigar string
	uint32_t *cigar = bam_get_cigar(read);
	cigarShift shift = ReadUtils::isUnmapped(read, hdr) ? cigarShift{nullptr, 0, 0} : hardClipCigar(cigar,
	                                                                                                read->core.n_cigar,
	                                                                                                start, stop);

	// the cigar may force a shift left or right (or both) in case we are left with insertions
	// starting or ending the read after applying the hard clip on start/stop.
	int newLength = read->core.l_qseq - (stop - start + 1) - shift.start - shift.end;

	// If the new read is going to be empty, return an empty read now. This avoids initializing the new
	// read with invalid values below in certain cases (such as a negative alignment start).
	// See https://github.com/broadinstitute/gatk/issues/3466
	if (newLength == 0)
		return nullptr;     //---to use nullptr or an empty bam1_t?

	char *newBases = new char[newLength + 1];
	char *newQuals = new char[newLength + 1];
	int copyStart = (start == 0) ? stop + 1 + shift.start : shift.start;

	// bam_set1() requires sequence representation of seq, so bases need to be decodedfrom 4 bit to ASCII code
	for (int i = 0; i < newLength; ++i) {
		uint8_t base = bam_seqi(bam_get_seq(read), copyStart + i);
		newBases[i] = ReadUtils::decodeBase(base);
	}
	memcpy(newQuals, bam_get_qual(read) + copyStart, newLength);
	newBases[newLength] = '\0';
	newQuals[newLength] = '\0';

	hts_pos_t newPosition = read->core.pos;

	bam1_t *hardClippedRead = bam_init1();
	if (start == 0 && !ReadUtils::isUnmapped(read, hdr)) {
		newPosition =
				read->core.pos + ReadUtils::calculateAlignmentStartShift(read->core.n_cigar, cigar, stop - start + 1);
	}
	//TODO: bam_set() can be omitted // there is a bug in bam_set1 and setClipBases owing to cigar calculated before this method
	bam_set1(hardClippedRead, read->core.l_qname, bam_get_qname(read), read->core.flag, read->core.tid, newPosition,
	         read->core.qual, shift.cigar->size(), shift.cigar->data(),
	         read->core.mtid, read->core.mpos, read->core.isize, newLength, newBases, newQuals, bam_get_l_aux(read));

	/*
		if(setClipBases(read, newPosition, shift.cigar->size(), shift.cigar->data(), newLength, newBases, newQuals) < 0)
		{
			throw "error in applyHardClipBases method";
		}
	*/
	/*
	uint8_t * qualities = bam_get_qual(hardClippedRead);
	for(int i=0; i<hardClippedRead->core.l_qseq; i++)
		cout << (int)qualities[i] << " ";
	cout << endl;
	*/
//    bam_destroy1(read);
//    read = hardClippedRead;

	delete shift.cigar;
	delete[] newBases;
	delete[] newQuals;
	return hardClippedRead;
}

bam1_t *ReadClipper::hardClipRead(bam1_t *read, int start, int stop, sam_hdr_t *hdr) {
	bam1_t *clippedRead = nullptr;
	int readLength = read->core.l_qseq;
	if (start < readLength) {
		if (stop >= readLength) {
			stop = readLength - 1;
		}
		clippedRead = applyHardClipBases(read, start, stop, hdr);
	}

	return clippedRead;
}
