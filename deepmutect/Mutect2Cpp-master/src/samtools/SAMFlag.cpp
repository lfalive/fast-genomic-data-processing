//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMFlag.h"

SAMFlag::SAMFlag(int flag, std::string& description) : flag(flag), description(description){

}

SAMFlag &SAMFlag::READ_PAIRED() {
    std::string tmp("Template having multiple segments in sequencing");
    static SAMFlag ret(1, tmp);
    return ret;
}

SAMFlag &SAMFlag::READ_UNMAPPED() {
    std::string tmp("Segment unmapped");
    static SAMFlag ret(4, tmp);
    return ret;
}

SAMFlag &SAMFlag::MATE_UNMAPPED() {
    std::string tmp("Next segment in the template unmapped");
    static SAMFlag ret(8, tmp);
    return ret;
}

SAMFlag &SAMFlag::PROPER_PAIR() {
    std::string tmp("Each segment properly aligned according to the aligner");
    static SAMFlag ret(2, tmp);
    return ret;
}

SAMFlag &SAMFlag::READ_REVERSE_STRAND() {
    std::string tmp("SEQ being reverse complemented");
    static SAMFlag ret(16, tmp);
    return ret;
}

SAMFlag &SAMFlag::MATE_REVERSE_STRAND() {
    std::string tmp("SEQ of the next segment in the template being reverse complemented");
    static SAMFlag ret(32, tmp);
    return ret;
}

SAMFlag &SAMFlag::FIRST_OF_PAIR() {
    std::string tmp("The first segment in the template");
    static SAMFlag ret(64, tmp);
    return ret;
}

SAMFlag &SAMFlag::SECOND_OF_PAIR() {
    std::string tmp("The last segment in the template");
    static SAMFlag ret(128, tmp);
    return ret;
}

SAMFlag &SAMFlag::SECONDARY_ALIGNMENT() {
    std::string tmp("Secondary alignment");
    static SAMFlag ret(256, tmp);
    return ret;
}

SAMFlag &SAMFlag::READ_FAILS_VENDOR_QUALITY_CHECK() {
    std::string tmp("Not passing quality controls");
    static SAMFlag ret(512, tmp);
    return ret;
}

SAMFlag &SAMFlag::DUPLICATE_READ() {
    std::string tmp("PCR or optical duplicate");
    static SAMFlag ret(1024, tmp);
    return ret;
}

SAMFlag &SAMFlag::SUPPLEMENTARY_ALIGNMENT() {
    std::string tmp("Supplementary alignment");
    static SAMFlag ret(2048, tmp);
    return ret;
}
