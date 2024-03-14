//
// Created by 梦想家xixi on 2021/10/13.
//

#ifndef MUTECT2CPP_MASTER_ACTIVITYPROFILESTATE_H
#define MUTECT2CPP_MASTER_ACTIVITYPROFILESTATE_H
#include "SimpleInterval.h"
#include "htslib/sam.h"
#include <iostream>

class SimpleInterval;

/**
 * The type of the value returned by {@link #getResultValue}
 */
enum Type {
    NONE,
    HIGH_QUALITY_SOFT_CLIPS
};

class ActivityProfileState {
private:
    SimpleInterval loc;
    double activeProb;

public:
    ActivityProfileState(SimpleInterval const &loc, double activeProb);

    ActivityProfileState(ActivityProfileState const &activityProfileState);

    ActivityProfileState(const std::string& refName, hts_pos_t pos, double activeProb);

    ActivityProfileState(int refName, hts_pos_t pos, double activeProb);

    virtual ~ActivityProfileState() = default;

    double isActiveProb() const;

    /**
     * Set the probability that this site is active.
     *
     * Probabilities should be between 0.0 and 1.0, however this is not currently enforced
     * because the {@link BandPassActivityProfile} can sometimes generate probabilities that
     * slightly exceed 1.0 when moving probability mass around. We intend to fix this by
     * capping at 1.0, but first we must evaluate the effects of capping on the HaplotypeCaller.
     *
     * @param activeProb probability (should be between 0.0 and 1.0) that the site is active
     */
    void setIsActiveProb(double activeProb);

    /**
    * The offset of state w.r.t. our current region's start location
    * @param regionStartLoc the start of the region, as a Locatable
    * @return the position of this profile relative to the start of this region
    */
    int getOffset(Locatable *regionStartLoc);

    SimpleInterval& getLoc();

    friend std::ostream & operator<<(std::ostream & os, ActivityProfileState& activityProfileState);
};


#endif //MUTECT2CPP_MASTER_ACTIVITYPROFILESTATE_H
