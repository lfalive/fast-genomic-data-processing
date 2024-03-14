//
// Created by lhh on 6/2/22.
//

#ifndef MUTECT2CPP_MASTER_VARIANTANNOTATION_H
#define MUTECT2CPP_MASTER_VARIANTANNOTATION_H

#include <string>
#include <vector>

/**
 * Superclass of all variant annotations.
 */
class VariantAnnotation {
public:
    virtual std::vector<std::string> getKeyNames() = 0;

    virtual std::string toString() = 0;
};


#endif //MUTECT2CPP_MASTER_VARIANTANNOTATION_H
