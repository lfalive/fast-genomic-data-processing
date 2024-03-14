//
// Created by 梦想家xixi on 2021/11/8.
//

#include <cstring>
#include <memory>
#include "Allele.h"
#include <stdexcept>
#include "StringUtils.h"

const std::shared_ptr<uint8_t[]> Allele::EMPTY_ALLELE_BASES(new uint8_t[0]);
const std::string Allele::NO_CALL_STRING = ".";
const std::string Allele::SPAN_DEL_STRING = "*";
const std::string Allele::NON_REF_STRING = "<NON_REF>";
const std::string Allele::UNSPECIFIED_ALTERNATE_ALLELE_STRING = "<*>";
std::shared_ptr<Allele> Allele::REF_A = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'A'}), 1,
                                                                 true);
std::shared_ptr<Allele> Allele::ALT_A = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'A'}), 1,
                                                                 false);
std::shared_ptr<Allele> Allele::REF_C = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'C'}), 1,
                                                                 true);
std::shared_ptr<Allele> Allele::ALT_C = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'C'}), 1,
                                                                 false);
std::shared_ptr<Allele> Allele::REF_G = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'G'}), 1,
                                                                 true);
std::shared_ptr<Allele> Allele::ALT_G = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'G'}), 1,
                                                                 false);
std::shared_ptr<Allele> Allele::REF_T = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'T'}), 1,
                                                                 true);
std::shared_ptr<Allele> Allele::ALT_T = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'T'}), 1,
                                                                 false);
std::shared_ptr<Allele> Allele::REF_N = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'N'}), 1,
                                                                 true);
std::shared_ptr<Allele> Allele::ALT_N = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'N'}), 1,
                                                                 false);
std::shared_ptr<Allele> Allele::SPAN_DEL = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'*'}), 1,
                                                                    false);
std::shared_ptr<Allele> Allele::NO_CALL = std::make_shared<Allele>(std::shared_ptr<uint8_t[]>(new uint8_t[1]{'.'}), 1,
                                                                   false);
std::shared_ptr<Allele> Allele::NON_REF_ALLELE = std::make_shared<Allele>(
		std::shared_ptr<uint8_t[]>(new uint8_t[9]{'<', 'N', 'O', 'N', '_', 'R', 'E', 'F', '>'}), 9, false);
std::shared_ptr<Allele> Allele::UNSPECIFIED_ALTERNATE_ALLELE = std::make_shared<Allele>(
		std::shared_ptr<uint8_t[]>(new uint8_t[3]{'<', '*', '>'}), 3, false);

Allele::Allele(std::shared_ptr<uint8_t[]> bases, int length, bool isRef) : isRef(false), bases(nullptr), length(0) {
	this->isRef = false;
	this->isNoCall = false;
	this->isSymbolic = false;
	this->bases = nullptr;
	if (wouldBeNullAllele(bases, length)) {
		throw std::invalid_argument("Null alleles are not supported");
	} else if (wouldBeNoCallAllele(bases, length)) {
		this->bases = EMPTY_ALLELE_BASES;
		this->length = length;
		this->isNoCall = true;
		if (isRef) {
			throw std::invalid_argument("Cannot tag a NoCall allele as the reference allele");
		}
	} else {
		if (wouldBeSymbolicAllele(bases, length)) {
			this->isSymbolic = true;
			if (isRef) {
				throw std::invalid_argument("Cannot tag a symbolic allele as the reference allele");
			}
		} else {
			StringUtils::toUpperCase(bases, length);
		}

		this->isRef = isRef;
		this->bases = bases;
		this->length = length;
		if (!acceptableAlleleBases(bases, length, isRef)) {
			throw std::invalid_argument("Unexpected base in allele bases");
		}
	}
}

bool Allele::wouldBeNullAllele(const std::shared_ptr<uint8_t[]> &bases, int length) {
	return length == 1 && bases.get()[0] == 45 || length == 0;
}

bool Allele::wouldBeNoCallAllele(const std::shared_ptr<uint8_t[]> &bases, int length) {
	return length == 1 && bases.get()[0] == 46;
}

bool Allele::wouldBeSymbolicAllele(const std::shared_ptr<uint8_t[]> &bases, int length) {
	if (length <= 1) {
		return false;
	} else {
		return bases.get()[0] == 60 || bases.get()[length - 1] == 62 || wouldBeBreakpoint(bases, length) ||
		       wouldBeSingleBreakend(bases, length);
	}
}

bool Allele::wouldBeBreakpoint(const std::shared_ptr<uint8_t[]> &bases, int length) {
	if (length <= 1) {
		return false;
	} else {
		for (int i = 0; i < length; ++i) {
			uint8_t base = bases.get()[i];
			if (base == 93 || base == 91) {
				return true;
			}
		}
		return false;
	}
}

bool Allele::wouldBeSingleBreakend(const std::shared_ptr<uint8_t[]> &bases, int length) {
	if (length <= 1) {
		return false;
	} else {
		return bases.get()[0] == 46 || bases.get()[length - 1] == 46;
	}
}

bool Allele::acceptableAlleleBases(const std::shared_ptr<uint8_t[]> &bases, int length, bool isReferenceAllele) {
	if (wouldBeNullAllele(bases, length)) {
		return false;
	} else if (!wouldBeNoCallAllele(bases, length) && !wouldBeSymbolicAllele(bases, length)) {
		if (wouldBeStarAllele(bases, length)) {
			return !isReferenceAllele;
		} else {
			std::shared_ptr<uint8_t[]> var2 = bases;
			int var3 = length;
			int var4 = 0;

			while (var4 < var3) {
				uint8_t base = var2.get()[var4];
				switch (base) {
					case 65:
					case 67:
					case 71:
					case 78:
					case 84:
					case 97:
					case 99:
					case 103:
					case 110:
					case 116:
						++var4;
						break;
					default:
						return false;
				}
			}

			return true;
		}
	} else {
		return true;
	}
}

bool Allele::wouldBeStarAllele(const std::shared_ptr<uint8_t[]> &bases, int length) {
	return length == 1 && bases.get()[0] == 42;
}

Allele::Allele(Allele &allele, bool ignoreRefState) : bases(allele.bases), isNoCall(allele.isNoCall),
                                                      isSymbolic(allele.isSymbolic) {
	this->isRef = !ignoreRefState && allele.isRef;
}

std::shared_ptr<Allele> Allele::create(std::shared_ptr<uint8_t[]> bases, int length, bool isRef) {
	if (bases == nullptr) {
		throw std::invalid_argument(
				"create: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");
	} else if (length == 1) {
		switch (bases.get()[0]) {
			case 42:
				if (isRef) {
					throw std::invalid_argument("Cannot tag a spanning deletions allele as the reference allele");
				}
				return SPAN_DEL;

			case 46:
				if (isRef) {
					throw std::invalid_argument("Cannot tag a NoCall allele as the reference allele");
				}
				return NO_CALL;

			case 65:
			case 97:
				return isRef ? REF_A : ALT_A;
			case 67:
			case 99:
				return isRef ? REF_C : ALT_C;
			case 71:
			case 103:
				return isRef ? REF_G : ALT_G;
			case 78:
			case 110:
				return isRef ? REF_N : ALT_N;
			case 84:
			case 116:
				return isRef ? REF_T : ALT_T;
			default:
				throw std::invalid_argument("Illegal base seen in the allele");
		}
	} else {
		return std::make_shared<Allele>(bases, length, isRef);
	}
}

std::shared_ptr<Allele> Allele::create(uint8_t base, bool isRef) {
	std::shared_ptr<uint8_t[]> bases(new uint8_t[1]{base});
	std::shared_ptr<Allele> res = create(bases, 1, isRef);
	return res;
}

std::shared_ptr<Allele> Allele::create(uint8_t base) {
	return create(base, false);
}

std::shared_ptr<Allele>
Allele::extend(const std::shared_ptr<Allele> &left, const std::shared_ptr<uint8_t[]> &right, int length) {
	if (left->isSymbolic) {
		throw std::invalid_argument("Cannot extend a symbolic allele");
	} else {
		std::shared_ptr<uint8_t[]> bases(new uint8_t[left->length + length]);
		memcpy(bases.get(), left->bases.get(), left->length);
		memcpy(bases.get() + left->length, right.get(), length);
		return create(bases, left->length + length, left->isRef);
	}
}

bool Allele::operator<(const Allele &other) const {
	if (isRef < other.getIsReference())
		return true;
	if (isRef > other.getIsReference())
		return false;
	if (isNoCall < other.getIsNoCall())
		return true;
	if (isNoCall < other.getIsNoCall())
		return false;
	return memcmp(bases.get(), other.bases.get(), length) <= 0;
}

int Allele::getLength() const {
	return isSymbolic ? 0 : length;
}

bool Allele::operator==(const Allele &other) const {
	if (this == &other)
		return true;
	if (this->isRef != other.getIsReference() || this->isNoCall != other.getIsNoCall())
		return false;
	if (this->getLength() != other.getLength())
		return false;
	if (this->bases == other.getBases())
		return true;
	return memcmp(bases.get(), other.bases.get(), length) == 0;
}

std::string Allele::getBaseString() {
	char tmp[length + 1];
	memcpy(tmp, bases.get(), length);
	tmp[length] = '\0';
	std::string ret(tmp, length);
	return ret;
}

size_t Allele::hashcode() {
	std::hash<std::string> hash_string;
	return hash_string(getBaseString());
}

bool Allele::equals(Allele &other, bool ignoreRefState) {
	if (this == &other)
		return true;
	if (isRef != other.isRef && !ignoreRefState)
		return false;
	if (isNoCall != other.isNoCall)
		return false;
	if (bases == other.bases)
		return true;
	if (length != other.length)
		return false;
	return memcmp(bases.get(), other.bases.get(), length) == 0;
}

bool Allele::getIsNonRefAllele() const{
    return this == NON_REF_ALLELE.get() || this == UNSPECIFIED_ALTERNATE_ALLELE.get();
}
