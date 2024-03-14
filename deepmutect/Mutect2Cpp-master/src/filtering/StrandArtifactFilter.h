//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_STRANDARTIFACTFILTER_H
#define MUTECT2CPP_MASTER_STRANDARTIFACTFILTER_H

#include "Mutect2VariantFilter.h"
#include "clustering/EStep.h"
#include "ErrorProbabilities.h"

class StrandArtifactFilter : public Mutect2VariantFilter{
private:
    constexpr static const double INITIAL_STRAND_ARTIFACT_PRIOR = 0.001;
    double strandArtifactPrior;
    double alphaStrand;
    double betaStrand;
    double INITIAL_ALPHA_STRAND;
    double INITIAL_BETA_STRAND;
    constexpr static double ALPHA_SEQ = 1;
    constexpr static double BETA_SEQ_SNV = 1000;
    constexpr static double BETA_SEQ_SHORT_INDEL = 5000;
    constexpr static double BETA_SEQ_LONG_INDEL = 50000;
    static const int LONG_INDEL_SIZE = 3;
    static const int LONGEST_STRAND_ARTIFACT_INDEL_SIZE = 4;
    std::vector<EStep> eSteps;

    EStep strandArtifactProbability(double strandArtifactPrior, int forwardCount, int reverseCount, int forwardAltCount, int reverseAltCount, int indelSize);
    double artifactStrandLogLikelihood(int strandCount, int strandAltCount);
    static double artifactStrandLogLikelihood(int strandCount, int strandAltCount, double alpha, double beta);
    double nonArtifactStrandLogLikelihood(int strandCount, int strandAltCount, int indelSize);

public:
    StrandArtifactFilter();
    virtual std::vector<std::string> requiredAnnotations();
    virtual double calculateErrorProbability(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine, std::shared_ptr<ReferenceContext>);
    EStep calculateArtifactProbabilities(const std::shared_ptr<VariantContext> & vc, Mutect2FilteringEngine* filteringEngine);
    virtual void accumulateDataForLearning(const std::shared_ptr<VariantContext> & vc, ErrorProbabilities errorProbabilities, Mutect2FilteringEngine* filteringEngine);
    virtual void learnParameters();
    virtual void clearAccumulatedData();
    virtual ErrorType errorType();
    virtual std::string filterName();
    int filterIndex() override;
};


#endif //MUTECT2CPP_MASTER_STRANDARTIFACTFILTER_H
