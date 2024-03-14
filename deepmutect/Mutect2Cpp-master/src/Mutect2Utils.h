//
// Created by 梦想家xixi on 2021/10/12.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2UTILS_H
#define MUTECT2CPP_MASTER_MUTECT2UTILS_H

#include <string>
#include <cfloat>
#include <vector>
#include <set>
#include <memory>

static double POSITIVE_INFINITY = DBL_MAX;

class Mutect2Utils
{
public:
    static std::string replaceWith(std::string& str1, const std::string& str2, const std::string& str3);
    static bool overlaps(int start, int end, int start2, int end2);
    static bool encloses(int outerStart, int outerEnd, int innerStart, int innerEnd);
    static void validateArg(bool condition, const std::string& msg);
    static bool goodProbability(double result);
    static double logLikelihoodRatio(int nRef, const std::vector<uint8_t>& altQuals, int repeatFactor);
    static double logLikelihoodRatio(int refCount, int altCount, double errorProbability);
    static int lastIndexOf(const std::shared_ptr<uint8_t[]>& reference, int refLength, const std::shared_ptr<uint8_t[]>& query, int queryLength);
    static std::shared_ptr<uint8_t[]> copyOfRange(const std::shared_ptr<uint8_t[]>& original, int ,int from, int to, int & length);
    static int Int_compare(int x, int y);
    static uint8_t decodeBase(uint8_t i);
    template<class T>
    static bool isVectorEquals(std::vector<T> a1, std::vector<T> a2){
        if(a1.size() != a2.size())
            return false;
        for(int i = 0; i < a1.size(); i++) {
            if(!((*a1.at(i)) == (*a2.at(i))))
                return false;
        }
        return true;
    }
    template<class T>
    static bool isSetEquals(std::set<T> a1, std::set<T> a2){
        if(a1.size() != a2.size())
            return false;
        typename std::set<T>::iterator iter1 = a1.begin();
        typename std::set<T>::iterator iter2 = a2.begin();
        while(iter1 != a1.end()) {
            if(!((**iter1) == (**iter2)))
                return false;
            iter1++;
            iter2++;
        }
        return true;
    }
};

#endif //MUTECT2CPP_MASTER_MUTECT2UTILS_H