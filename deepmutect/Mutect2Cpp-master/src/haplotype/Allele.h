//
// Created by 梦想家xixi on 2021/11/8.
//

#ifndef MUTECT2CPP_MASTER_ALLELE_H
#define MUTECT2CPP_MASTER_ALLELE_H


#include <cstdint>
#include <string>
#include <memory>
#include <xxhash.hpp>

class Allele {
private:
    static const std::shared_ptr<uint8_t[]> EMPTY_ALLELE_BASES;
    static const char SINGLE_BREAKEND_INDICATOR = '.';
    static const char BREAKEND_EXTENDING_RIGHT = '[';
    static const char BREAKEND_EXTENDING_LEFT = ']';
    static const char SYMBOLIC_ALLELE_START = '<';
    static const char SYMBOLIC_ALLELE_END = '>';
    bool isRef;
    bool isNoCall;
    bool isSymbolic;
    std::shared_ptr<uint8_t[]> bases;
    int length;

public:
    static std::shared_ptr<Allele> REF_A;
    static std::shared_ptr<Allele> ALT_A;
    static std::shared_ptr<Allele> REF_C;
    static std::shared_ptr<Allele> ALT_C;
    static std::shared_ptr<Allele> REF_G;
    static std::shared_ptr<Allele> ALT_G;
    static std::shared_ptr<Allele> REF_T;
    static std::shared_ptr<Allele> ALT_T;
    static std::shared_ptr<Allele> REF_N;
    static std::shared_ptr<Allele> ALT_N;
    static std::shared_ptr<Allele> SPAN_DEL;
    static std::shared_ptr<Allele> NO_CALL;
    static std::shared_ptr<Allele> NON_REF_ALLELE;
    static std::shared_ptr<Allele> UNSPECIFIED_ALTERNATE_ALLELE;
    static Allele SV_SIMPLE_DEL;
    static Allele SV_SIMPLE_INS;
    static Allele SV_SIMPLE_INV;
    static Allele SV_SIMPLE_CNV;
    static Allele SV_SIMPLE_DUP;

    static const std::string NO_CALL_STRING;
    static const std::string SPAN_DEL_STRING;
    static const std::string NON_REF_STRING;
    static const std::string UNSPECIFIED_ALTERNATE_ALLELE_STRING;
    static bool wouldBeNullAllele(const std::shared_ptr<uint8_t[]> & bases, int length);
    static bool wouldBeNoCallAllele(const std::shared_ptr<uint8_t[]> & bases, int length);
    static bool wouldBeSymbolicAllele(const std::shared_ptr<uint8_t[]> & bases, int length);
    static bool wouldBeBreakpoint(const std::shared_ptr<uint8_t[]> & bases, int length);
    static bool wouldBeSingleBreakend(const std::shared_ptr<uint8_t[]> & bases, int length);
    static bool wouldBeStarAllele(const std::shared_ptr<uint8_t[]> & bases, int length);
    static bool acceptableAlleleBases(const std::shared_ptr<uint8_t[]> & bases, int length, bool isReferenceAllele);
    static std::shared_ptr<Allele> create(std::shared_ptr<uint8_t[]> bases, int length, bool isRef);
    static std::shared_ptr<Allele> create(uint8_t base, bool isRef);
    static std::shared_ptr<Allele> create(uint8_t base);
    static std::shared_ptr<Allele> extend(const std::shared_ptr<Allele> & left, const std::shared_ptr<uint8_t[]> & right, int length);
    bool getIsNoCall() const  {return isNoCall;}
    bool getIsCalled() const {return !isNoCall;}
    bool getIsReference() const {return isRef;}
    bool getIsNonReference() const {return !isRef;}
    bool getIsNonRefAllele() const;
    bool getIsSymbolic() const {return isSymbolic;}
    bool getIsBreakpoint() const {return wouldBeBreakpoint(bases, length);}
    bool getIsSingleBreakend() const {return wouldBeSingleBreakend(bases, length);}
    std::shared_ptr<uint8_t[]> getBases() const {return bases;}
    bool equals(Allele & other, bool ignoreRefState);
    bool operator<(const Allele & other) const;
    bool operator==(const Allele & other) const;
    int getLength() const;
    int getBasesLength() const {return length;}
    std::string getBaseString();
    size_t hashcode();
    Allele(std::shared_ptr<uint8_t[]> bases, int length, bool isRef);
    Allele(Allele &  allele, bool ignoreRefState);

protected:

};

struct hash_Allele {
    size_t operator()(const std::shared_ptr<Allele> &allele) const {
        return xxh::xxhash3<64>(allele->getBases().get(), allele->getBasesLength());
    }
};

struct equal_Allele {
    bool operator()(const std::shared_ptr<Allele> &left, const std::shared_ptr<Allele> &right) const {
        if (left->getLength() != right->getLength())
            return false;
		return memcmp(left->getBases().get(), right->getBases().get(), left->getLength()) == 0;
    }
};



#endif //MUTECT2CPP_MASTER_ALLELE_H
