//
// Created by 梦想家xixi on 2021/10/11.
//

#include "SimpleInterval.h"
#include "boost/utility.hpp"
#include <iostream>

SimpleInterval::SimpleInterval(const std::string &contig, int start, int end) : contig(ContigMap::getContigInt(contig)),
                                                                                start(start), end(end) {
	if (BOOST_UNLIKELY(this->contig == -1) || start < 0 || start > end)
		throw std::invalid_argument("SimpleInterval::validatePositions failed.");
}

SimpleInterval::SimpleInterval(std::string &&contig, int start, int end) : contig(ContigMap::getContigInt(contig)),
                                                                           start(start), end(end) {}

SimpleInterval::SimpleInterval(int contig, int start, int end) : contig(contig), start(start), end(end) {
	if (BOOST_UNLIKELY(contig == -1) || start < 0 || start > end)
		throw std::invalid_argument("SimpleInterval::validatePositions failed.");
}

SimpleInterval::SimpleInterval(SimpleInterval const &simpleInterval) : contig(simpleInterval.contig),
                                                                       start(simpleInterval.start),
                                                                       end(simpleInterval.end) {}

SimpleInterval::SimpleInterval(const std::shared_ptr<Locatable> &pLocatable)
		: contig(ContigMap::getContigInt(pLocatable->getContig())), start(pLocatable->getStart()),
		  end(pLocatable->getEnd()) {}

SimpleInterval::SimpleInterval() : contig(-1), start(0), end(0) {}

SimpleInterval::SimpleInterval(std::string &str) {
	if (BOOST_UNLIKELY(str.empty()))
		throw std::invalid_argument("Null object is not allowed here.");

	std::string m_contig;
	int m_start;
	int m_end;

	const int colonIndex = (int) str.find_last_of(CONTIG_SEPARATOR);
	if (colonIndex == -1) {
		m_contig = str;
		m_start = 1;
		m_end = INT32_MAX;
	} else {
		m_contig = str.substr(0, colonIndex);
		const int dashIndex = (int) str.find(START_END_SEPARATOR, colonIndex);
		if (dashIndex == -1) {
			if (str.find_last_of(END_OF_CONTIG) == str.size() - END_OF_CONTIG.size()) {
				std::string pos = str.substr(colonIndex + 1, str.size() - colonIndex - 2);
				m_start = parsePosition(pos);
				m_end = INT32_MAX;
			} else {
				std::string pos = str.substr(colonIndex + 1, str.size() - colonIndex - 1);
				m_start = parsePosition(pos);
				m_end = m_start;
			}
		} else {
			std::string pos = str.substr(colonIndex + 1, dashIndex - colonIndex - 1);
			m_start = parsePosition(pos);
			pos = str.substr(dashIndex + 1, str.size() - dashIndex - 1);
			m_end = parsePosition(pos);
		}
	}

	int contigInt = ContigMap::getContigInt(m_contig);
	if (BOOST_UNLIKELY(contigInt == -1) || start < 0 || start > end)
		throw std::invalid_argument("SimpleInterval::validatePositions failed.");

	contig = contigInt;
	start = m_start;
	end = m_end;
}

void SimpleInterval::clearContig() {
	this->contig = -1;
}

int SimpleInterval::parsePosition(std::string pos) {
	int postion;
	pos = Mutect2Utils::replaceWith(pos, ",", "");
	std::stringstream ss;
	ss << pos;
	ss >> postion;
	if (ss.eof() && !ss.fail())
		return postion;

	throw std::invalid_argument("Contig input error.");
}

bool SimpleInterval::operator==(const SimpleInterval &interval) const {
	if (contig == interval.contig && start == interval.start && end == interval.end)
		return true;
	return false;
}

bool SimpleInterval::overlapsWithMargin(const std::shared_ptr<SimpleInterval> &other, const int margin) const {
	if (BOOST_UNLIKELY(margin < 0))
		throw std::invalid_argument("Given margin is negative.");

	if (BOOST_UNLIKELY(other == nullptr) || BOOST_UNLIKELY(other->getContigInt() == -1))
		return false;

	return this->contig == other->getContigInt() && this->start <= other->getEnd() + margin &&
	       other->getStart() - margin <= this->end;
}

bool SimpleInterval::overlaps(const std::shared_ptr<SimpleInterval> &other) const {
	return overlapsWithMargin(other, 0);
}

std::shared_ptr<SimpleInterval> SimpleInterval::intersect(const std::shared_ptr<SimpleInterval> &other) const {
	if (BOOST_LIKELY(overlaps(other)))
		return std::make_shared<SimpleInterval>(contig, std::max(start, other->getStart()),
		                                        std::min(end, other->getEnd()));

	throw std::invalid_argument("SimpleInterval::intersect(): The two intervals need to overlap.");
}

std::shared_ptr<SimpleInterval> SimpleInterval::mergeWithContiguous(const std::shared_ptr<SimpleInterval> &other) {
	if (BOOST_UNLIKELY(other == nullptr))
		throw std::invalid_argument("Null object is not allowed here.");

	if (BOOST_LIKELY(contiguous(other.get())))
		return std::make_shared<SimpleInterval>(contig, std::min(start, other->getStart()),
		                                        std::max(end, other->getEnd()));

	throw std::invalid_argument("The two intervals need to be contiguous.");
}

bool SimpleInterval::contiguous(SimpleInterval *other) const {
	return contig == other->getContigInt() && start <= other->getEnd() + 1 && other->getStart() <= end + 1;
}

std::shared_ptr<SimpleInterval> SimpleInterval::spanWith(const std::shared_ptr<SimpleInterval> &other) const {
	if (other == nullptr)
		throw std::invalid_argument("Null object is not allowed here.");
	if (contig != other->getContigInt())
		throw std::invalid_argument("Cannot get span for intervals on different contigs.");

	return std::make_shared<SimpleInterval>(contig, std::min(start, other->getStart()),
	                                        std::max(end, other->getEnd()));
}

std::shared_ptr<SimpleInterval> SimpleInterval::expandWithinContig(const int padding, const int contigLength) const {
	if (padding < 0)
		throw std::invalid_argument("Padding must be >= 0.");

	return IntervalUtils::trimIntervalToContig(contig, start - padding, end + padding,
	                                           contigLength);
}

std::ostream &operator<<(std::ostream &os, const SimpleInterval &simpleInterval) {
	os << "contig:" << simpleInterval.contig << "  start:" << simpleInterval.start << "   end:" << simpleInterval.end
	   << std::endl;
	return os;
}

std::shared_ptr<SimpleInterval>
SimpleInterval::expandWithinContig(int padding, SAMSequenceDictionary *sequenceDictionary) const {
	if (BOOST_UNLIKELY(sequenceDictionary == nullptr))
		throw std::invalid_argument("null is not allowed there");

	SAMSequenceRecord &contigRecord = sequenceDictionary->getSequences()[contig];
	return expandWithinContig(padding, contigRecord.getSequenceLength());
}

void SimpleInterval::printInfo() const {
	std::cout << getContig() << " " << getStart() + 1 << " " << getEnd() + 1 << std::endl;
}

