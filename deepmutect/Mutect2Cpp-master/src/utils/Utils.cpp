//
// Created by lhh on 4/26/22.
//

#include <cstring>
#include <sstream>
#include "Utils.h"

bool Utils::equalRange(uint8_t *left, int leftOffset, uint8_t *right, int rightOffset, int length)
{
    for(int i=0; i<length; i++)
    {
        if(left[leftOffset + i] != right[rightOffset + i])
            return false;
    }
    return true;
}

std::shared_ptr<char[]> Utils::dupBytes(char b, int nCopies)
{
    std::shared_ptr<char[]> bytes(new char[nCopies]);
    memset(bytes.get(), (int)b, nCopies);
    return bytes;
}

std::string Utils::join(std::string&& separator, std::set<std::string>& objects) {
    if(objects.empty())
        return "";

    auto iter = objects.begin();
    if(objects.size() == 1) // fast path for singleton collections
        return *iter;
    else{   // full path for 2+ collection that actually need a join
        std::stringstream ret(*iter);
        while(iter != objects.end())
        {
            iter++;
            ret << separator;
            ret << *iter;
        }
        return ret.str();
    }
}