//
// Created by 梦想家xixi on 2021/12/20.
//

#ifndef MUTECT2CPP_MASTER_SAMFLAG_H
#define MUTECT2CPP_MASTER_SAMFLAG_H

#include <string>

class SAMFlag {
private:
    SAMFlag(int flag, std::string& description);

public:
    const int flag;
    std::string  description;
    static SAMFlag & READ_PAIRED();
    static SAMFlag & PROPER_PAIR();
    static SAMFlag & READ_UNMAPPED();
    static SAMFlag & MATE_UNMAPPED();
    static SAMFlag & READ_REVERSE_STRAND();
    static SAMFlag & MATE_REVERSE_STRAND();
    static SAMFlag & FIRST_OF_PAIR();
    static SAMFlag & SECOND_OF_PAIR();
    static SAMFlag &SECONDARY_ALIGNMENT();
    static SAMFlag & READ_FAILS_VENDOR_QUALITY_CHECK();
    static SAMFlag & DUPLICATE_READ();
    static SAMFlag &SUPPLEMENTARY_ALIGNMENT();
};


#endif //MUTECT2CPP_MASTER_SAMFLAG_H
