//
// The class used to represent attribute value
// Created by lhh on 6/2/22.
//

#ifndef MUTECT2CPP_MASTER_ATTRIBUTEVALUE_H
#define MUTECT2CPP_MASTER_ATTRIBUTEVALUE_H

#include <vector>
#include <string>

enum Attribute_Type{
    INT,
    DOUBLE,
    BOOL,
    STRING,
    VECTOR_DOUBLE,
    VECTOR_INT,
    NULL_       // represent empty Attribute value
};


class AttributeValue {
private:
    void* value;
    Attribute_Type type;

public:


    AttributeValue(void* value, Attribute_Type type);

    AttributeValue(int val);

    AttributeValue(double val);

    AttributeValue(bool val);

    AttributeValue(std::string& val);

    AttributeValue(std::vector<double>& val);

    AttributeValue(std::vector<int>& val);

    AttributeValue(const AttributeValue& other);

    ~AttributeValue();

    int getAttributeAsInt() const;

    double getAttributeAsDouble();

    std::vector<double> getAttributeAsDoubleVector() const;

    std::vector<int> getAttributeAsIntVector() const;

    std::string getAttributeAsString() const;

    static AttributeValue empty_value();

    bool getAttributeAsBool() const;
};


#endif //MUTECT2CPP_MASTER_ATTRIBUTEVALUE_H
