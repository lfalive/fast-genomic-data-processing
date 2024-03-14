//
// Created by lhh on 10/25/21.
//

#include "BandPassActivityProfile.h"
#include "MathUtils.h"

BandPassActivityProfile::BandPassActivityProfile(int maxProbPropagationDistance, double activeProbThreshold,
                                                 int maxFilterSize, double sigma, bool adaptiveFilterSize,
                                                 SAMFileHeader *header) : ActivityProfile(maxProbPropagationDistance, activeProbThreshold, header), sigma(sigma)
{
    std::vector<double> * fullKernel = makeKernel(maxFilterSize, sigma);
    filterSize = adaptiveFilterSize ? determineFilterSize(fullKernel, MIN_PROB_TO_KEEP_IN_FILTER) : maxFilterSize;
    gaussianKernel = makeKernel(filterSize, sigma);
    negativeState = new std::vector<std::shared_ptr<ActivityProfileState>>();
    negativeState->reserve(1);
    activateState = new std::vector<std::shared_ptr<ActivityProfileState>>();
    activateState->reserve(filterSize*2 + 1);
    delete fullKernel;
}

std::vector<double> * BandPassActivityProfile::makeKernel(int filterSize, double sigma)
{
    int bandSize = (filterSize << 1) + 1;
    auto * kernel = new std::vector<double>;
    kernel->reserve(bandSize);
    for(int i=0; i<bandSize; i++)
    {
        kernel->push_back(MathUtils::normalDistribution(filterSize, sigma, i));   //---to be tested
    }
    return MathUtils::normalizeSumToZero(kernel);
}

BandPassActivityProfile::~BandPassActivityProfile()
{
    delete gaussianKernel;
    delete negativeState;
    delete activateState;
}

int BandPassActivityProfile::determineFilterSize(std::vector<double> * kernel, double minProbToKeepInFilter)
{
    int middle = (kernel->size() - 1) >> 1;
    int filterEnd = middle;
    while (filterEnd > 0){
        if(kernel->operator[](filterEnd - 1) < minProbToKeepInFilter)
            break;
        filterEnd--;
    }
    return middle - filterEnd;
}

vector<std::shared_ptr<ActivityProfileState>> * BandPassActivityProfile::processState(const std::shared_ptr<ActivityProfileState> & justAddedState)
{
    double activeProb = justAddedState->isActiveProb();

    if(activeProb > 0.0)
    {
        for(int i = -filterSize; i <= filterSize; i++)
        {
            optional<SimpleInterval> loc = getLocForOffset(justAddedState->getLoc(), i);
            if(loc.has_value())
            {
                double newProb = activeProb * gaussianKernel->operator[](i + filterSize);
                const std::shared_ptr<ActivityProfileState> profile = std::make_shared<ActivityProfileState>(loc.value(), newProb);
                activateState->emplace_back(profile);

            }
        }
        return activateState;
    } else{
        negativeState->emplace_back(justAddedState);
        return negativeState;
    }
}

int BandPassActivityProfile::getMaxProbPropagationDistance()
{
    return maxProbPropagationDistance + filterSize;
}


