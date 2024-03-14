//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMBinaryTagAndValue.h"
#include "SAMUtils.h"
#include <stdexcept>
#include <iostream>

SAMBinaryTagAndValue::SAMBinaryTagAndValue(short tag, void *value, Void_Type voidType, int length) : tag(tag), value(value), type(voidType){
    if(value == nullptr) {
        throw std::invalid_argument("SAMBinaryTagAndValue value may not be null");
    } else if (!isAllowedAttributeValue(value, voidType)) {
        throw std::invalid_argument("Attribute type not supported.");
    } else {
        if((voidType != Uint8_t_Array_Type) && (voidType != Short_Array_Type) && (voidType != Int_Array_Type) && (voidType != Float_Array_Type))
            this->length = -1;
        else
            this->length = length;
    }
}

bool SAMBinaryTagAndValue::isAllowedAttributeValue(void* value, Void_Type voidType) {
    if((voidType != Uint8_Type) && (voidType != Short_Type) && (voidType != Integer_Type) && (voidType != String_Type) && (voidType != Character_Type)
    && (voidType != Float_Type) && (voidType != Uint8_t_Array_Type) && (voidType != Short_Array_Type) && (voidType != Int_Array_Type) && (voidType != Float_Array_Type)) {
        if((voidType != Long_Type)) {
            return false;
        } else {
            return SAMUtils::isValidUnsignedIntegerAttribute(*(long*)value) || *(long*)value >= -2147483648L && *(long*)value <= 2147483647L;
        }
    } else {
        return true;
    }
}

std::shared_ptr<SAMBinaryTagAndValue> SAMBinaryTagAndValue::remove(std::shared_ptr<SAMBinaryTagAndValue> root,short tag) {
    if(root->tag == tag){
        std::shared_ptr<SAMBinaryTagAndValue> ret = root->next;
        return ret;
    }
    std::shared_ptr<SAMBinaryTagAndValue> iter = root;
    while(iter->next != nullptr) {
        if(tag == iter->next->tag) {
            std::shared_ptr<SAMBinaryTagAndValue> tmp = iter->next;
            iter->next = tmp->next;
            break;
        } else {
            iter = iter->next;
        }
    }
    return root;
}

std::shared_ptr<SAMBinaryTagAndValue> SAMBinaryTagAndValue::insert(std::shared_ptr<SAMBinaryTagAndValue> root, std::shared_ptr<SAMBinaryTagAndValue> attr) {
    if(attr == nullptr)
        return root;
    else if (attr->next != nullptr) {
        throw std::invalid_argument("Can only insert single tag/value combinations.");
    }

    if (attr->tag < root->tag){
        // attr joins the list ahead of this element
        attr->next = root;
        return attr;
    } else if(attr->tag == root->tag){
        // attr replaces this in the list
        attr->next = root->next;
        return attr;
    } else if(root->next == nullptr) {
        // attr gets stuck on the end
        root->next = attr;
        return root;
    } else {
        // attr gets inserted somewhere in the tail
        root->next = root->next->insert(root->next, attr);
        return root;
    }
}

SAMBinaryTagAndValue* SAMBinaryTagAndValue::find(short tag) {
    if(this->tag == tag)
        return this;
    else if(this->tag > tag || this->next == nullptr)
        return nullptr;
    else
        return this->next->find(tag);
}

SAMBinaryTagAndValue::~SAMBinaryTagAndValue() {
    delete (std::string *)value;
}
