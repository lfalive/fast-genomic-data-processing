//
// Created by 梦想家xixi on 2021/12/20.
//

#ifndef MUTECT2CPP_MASTER_SAMBINARYTAGANDVALUE_H
#define MUTECT2CPP_MASTER_SAMBINARYTAGANDVALUE_H

#include <memory>

enum Void_Type{
    Uint8_Type,
    Long_Type,
    Short_Type,
    Integer_Type,
    String_Type,
    Character_Type,
    Float_Type,
    Uint8_t_Array_Type,
    Short_Array_Type,
    Int_Array_Type,
    Float_Array_Type,
    NuLL_Type
};

class SAMBinaryTagAndValue {
public:
    short tag;
    void* value;
    Void_Type type;
    int length;
    std::shared_ptr<SAMBinaryTagAndValue> next = nullptr;
    SAMBinaryTagAndValue(short tag, void* value, Void_Type voidType, int length);
    static std::shared_ptr<SAMBinaryTagAndValue> remove(std::shared_ptr<SAMBinaryTagAndValue> root, short tag);
    static std::shared_ptr<SAMBinaryTagAndValue> insert(std::shared_ptr<SAMBinaryTagAndValue> root, std::shared_ptr<SAMBinaryTagAndValue> attr);
    SAMBinaryTagAndValue* find(short tag);
    ~SAMBinaryTagAndValue();

protected:
    static bool isAllowedAttributeValue(void* value, Void_Type type);
};


#endif //MUTECT2CPP_MASTER_SAMBINARYTAGANDVALUE_H
