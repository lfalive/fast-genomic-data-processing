//
// Created by lhh on 10/25/21.
//

#ifndef MUTECT2CPP_MASTER_BANDPASSACTIVITYPROFILE_H
#define MUTECT2CPP_MASTER_BANDPASSACTIVITYPROFILE_H

#include "ActivityProfile.h"
using namespace std;

/**
 * A band pass filtering of the activity profile
 *
 * Applies a band pass filter with a Gaussian kernel to the input state probabilities to smooth
 * them out of an interval
 */
class BandPassActivityProfile : public ActivityProfile{
private:
    constexpr static double MIN_PROB_TO_KEEP_IN_FILTER = 1e-5;

    int filterSize;
    double sigma;
    vector<double> * gaussianKernel;

    vector<double> * makeKernel(int filterSize, double sigma);

    std::vector<std::shared_ptr<ActivityProfileState>> * activateState;
    std::vector<std::shared_ptr<ActivityProfileState>> * negativeState;

    int determineFilterSize(vector<double> * kernel, double minProbToKeepInFilter);

public:
    const static int MAX_FILTER_SIZE = 50;
    constexpr static double DEFAULT_SIGMA = 17.0;

    BandPassActivityProfile(int maxProbPropagationDistance, double activeProbThreshold, int maxFilterSize,
                                double sigma, bool adaptiveFilterSize, SAMFileHeader * header);

    ~BandPassActivityProfile();

    /**
     * Band pass the probabilities in the ActivityProfile, producing a new profile that's band pass filtered
     * @return a new double[] that's the band-pass filtered version of this profile
     */
    vector<std::shared_ptr<ActivityProfileState>> * processState(const std::shared_ptr<ActivityProfileState> & justAddedState) override;

    int getMaxProbPropagationDistance() override;
};


#endif //MUTECT2CPP_MASTER_BANDPASSACTIVITYPROFILE_H
