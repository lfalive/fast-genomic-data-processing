//
// Created by 梦想家xixi on 2021/12/27.
//

#include <stdexcept>
#include "HeaderRecordFactory.h"

AbstractSAMHeaderRecord *HeaderRecordFactory::createRecord(std::string &newId, AbstractSAMHeaderRecord *record) {
    throw std::invalid_argument("function not implement");
}
