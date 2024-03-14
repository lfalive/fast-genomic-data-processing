//
// Created by lhh on 6/2/22.
//

#include <iostream>
#include <cassert>
#include "AttributeValue.h"


AttributeValue::AttributeValue(void *value, Attribute_Type type): value(value), type(type) {

}

AttributeValue::AttributeValue(int val) {
    value = (void*)(new int(val));
    type = INT;
}

AttributeValue::AttributeValue(double val) {
    value = (void*)(new double(val));
    type = DOUBLE;
}

AttributeValue::AttributeValue(bool val) {
    value = (void *)(new bool(val));
    type = BOOL;
}

AttributeValue::AttributeValue(std::string& val) {
    value = (void *)(new std::string(val));
    type = STRING;
}

AttributeValue::AttributeValue(std::vector<double> &val) {
    value = (void*)(new std::vector<double>(val));
    type = VECTOR_DOUBLE;
}

AttributeValue::AttributeValue(std::vector<int> &val) {
    value = (void*)(new std::vector<int>(val));
    type = VECTOR_INT;
}

AttributeValue::AttributeValue(const AttributeValue &other) {
    this->type = other.type;
    switch (type) {
        case INT:
            value = (void*)(new int(*(int*)(other.value)));
            break;
        case DOUBLE:
            value = (void*)(new double(*(double*)(other.value)));
            break;
        case BOOL:
            value = (void*)(new bool(*(bool*)(other.value)));
            break;
        case STRING:
            value = (void*)(new std::string(*(std::string*)(other.value)));
            break;
        case VECTOR_DOUBLE:
            value = (void*)(new std::vector<double>(*(std::vector<double>*)(other.value)));
            break;
        case VECTOR_INT:
            value = (void*)(new std::vector<int>(*(std::vector<int>*)(other.value)));
            break;
    }
}

int AttributeValue::getAttributeAsInt() const {
    assert(type == INT);
    return *((int*)value);
}

AttributeValue::~AttributeValue(){
    if(value != nullptr && type != NULL_)
    {
        switch (type) {
            case INT:
                delete (int*)value;
                break;
            case DOUBLE:
                delete (double*)value;
                break;
            case BOOL:
                delete (bool*)value;
                break;
            case STRING:
                delete (std::string* )value;
                break;
            case VECTOR_DOUBLE:
                delete (std::vector<double>*)value;
                break;
            case VECTOR_INT:
                delete (std::vector<int>*)value;
                break;
        }
    }
}

AttributeValue AttributeValue::empty_value() {
    return AttributeValue(nullptr, NULL_);
}

double AttributeValue::getAttributeAsDouble() {
    assert(type == DOUBLE);
    return *((double*)value);
}

std::vector<double> AttributeValue::getAttributeAsDoubleVector() const {
    assert(type == VECTOR_DOUBLE);
    return *((std::vector<double>*)value);
}

std::vector<int> AttributeValue::getAttributeAsIntVector() const {
    assert(type == VECTOR_INT);
    return *((std::vector<int>*)value);
}

std::string AttributeValue::getAttributeAsString() const{
    assert(type == STRING);
    return *((std::string*)value);
}

bool AttributeValue::getAttributeAsBool() const {
    assert(type == BOOL);
    return *((bool*)value);
}
