//
// Created by 梦想家xixi on 2021/11/8.
//

#ifndef MUTECT2CPP_MASTER_STRINGUTILS_H
#define MUTECT2CPP_MASTER_STRINGUTILS_H


#include <cstdint>
#include <string>
#include <vector>
#include <memory>

class StringUtils {
public:
    static void toUpperCase(std::shared_ptr<uint8_t[]> & bytes, int length);
    static double parseDouble(const std::string & str);
    static int parseInt(const std::string & str);
    static std::string doubleToString(double d);
    static std::string join(std::string &a, std::vector<std::string> &lists);
};


#endif //MUTECT2CPP_MASTER_STRINGUTILS_H
