//
// Created by lhh on 11/13/21.
//

#ifndef MUTECT2CPP_MASTER_COMMONINFO_H
#define MUTECT2CPP_MASTER_COMMONINFO_H

#include <memory>
#include <string>
#include <set>
#include <map>
#include "AttributeValue.h"

/*enum Attribute_Type{
    INT,
    DOUBLE,
    BOOL,
    STRING,
    VECTOR_DOUBLE
};*/

/**
 * Common utility routines for VariantContext and Genotype
 */
class CommonInfo {
private:
    double log10PError = NO_LOG10_PERROR;
    std::string name;
    std::set<std::string> filters;
    std::shared_ptr<std::map<std::string, AttributeValue>> attributes;

    /* 1:int
     * 2:double
     * 3:bool
     * 4:std::string
     */
    //std::map<std::string, Attribute_Type> attributeTotypeMap;   // Maybe this parameter is useless

public:
    constexpr static double NO_LOG10_PERROR = 1.0;

    CommonInfo(std::string & name, double log10PError, std::set<std::string> * filters, std::shared_ptr<std::map<std::string, AttributeValue>> attributes);

    ~CommonInfo();

    void setLog10PError(double log10PError);

    bool hasAttribute(const std::string &key);

    int getAttributeAsInt(std::string &key, int defaultValue);

    std::vector<int> getAttributeAsIntVector(const std::string &key, std::vector<int> defaultValue);

    AttributeValue getAttribute(std::string &key);

    const std::map<std::string, AttributeValue> & getAttributes();

    std::shared_ptr<std::map<std::string, AttributeValue>> getAttributesAsPointer();

    std::set<std::string> * getFiltersMaybeNull();

    std::set<std::string>& getFilters();

    double getLog10PError() const;

    std::string & getName();

    bool isNotFiltered(){
        return !isFiltered();
    }

    bool isFiltered() {
        return !filters.empty();
    }

    bool filtersWereApplied();


};


#endif //MUTECT2CPP_MASTER_COMMONINFO_H
