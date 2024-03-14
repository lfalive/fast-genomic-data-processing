//
// Created by 梦想家xixi on 2021/11/8.
//

#include "StringUtils.h"
#include <sstream>
#include <vector>

void StringUtils::toUpperCase(std::shared_ptr<uint8_t[]> & bytes_, int length) {
    uint8_t * bytes = bytes_.get();
    for(int i = 0; i < length; ++i) {
        if (bytes[i] >= 97 && bytes[i] <= 122) {
            bytes[i] += -32;
        }
    }
}

double StringUtils::parseDouble(const std::string &str) {
    double ret;
    std::stringstream ss;
    ss << str;
    ss >> ret;
    if(ss.eof() && !ss.fail())
        return ret;
    else
        throw std::invalid_argument("input string error");
}

int StringUtils::parseInt(const std::string &str) {
    int ret;
    std::stringstream ss;
    ss << str;
    ss >> ret;
    if(ss.eof() && !ss.fail())
        return ret;
    else
        throw std::invalid_argument("input string error");
}

std::string StringUtils::doubleToString(double d) {
    std::string ret;
    std::stringstream ss;
    ss.precision(15);
    ss << d;
    ss >> ret;
    if(ss.eof() && !ss.fail())
        return ret;
    else
        throw std::invalid_argument("input string error");
}

std::string StringUtils::join(std::string &a, std::vector<std::string> &lists) {
    std::string ret;
    std::stringstream ss;
    ss << a;
    for(const std::string& str : lists) {
        ss << str;
    }
    ss >> ret;
    return ret;
}
