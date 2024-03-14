//
// Created by lhh on 4/26/22.
//

#ifndef MUTECT2CPP_MASTER_UTILS_H
#define MUTECT2CPP_MASTER_UTILS_H

#include <cstdint>
#include <memory>
#include <set>
#include <string>

class Utils {
public:
    static bool equalRange(uint8_t* left, int leftOffset, uint8_t* right, int rightOffset, int length);

    static std::shared_ptr<char[]> dupBytes(char b, int nCopies);

    /**
   * Returns a string of the form elt1.toString() [sep elt2.toString() ... sep elt.toString()] for a collection of
   * elti objects (note there's no actual space between sep and the elti elements).  Returns
   * "" if collection is empty.  If collection contains just elt, then returns elt.toString()
   *
   * @param separator the string to use to separate objects
   * @param objects a collection of objects.  the element order is defined by the iterator over objects
   * @param <T> the type of the objects
   * @return a non-null string
   */
    static std::string join(std::string&& separator, std::set<std::string>& objects);
};


#endif //MUTECT2CPP_MASTER_UTILS_H
