//
// Created by lhh on 5/28/22.
//

#ifndef MUTECT2CPP_MASTER_SOMATICLIKELIHOODSENGINE_H
#define MUTECT2CPP_MASTER_SOMATICLIKELIHOODSENGINE_H

#include <vector>
#include <memory>

class SomaticLikelihoodsEngine {
private:
    constexpr static double NEGLIGIBLE_RESPONSIBILITY = 1.0e-10;



    static double xLogx(double x);

    static double likelihoodsContribution(const std::shared_ptr<std::vector<double>>& logLikelihoodsForRead, const std::shared_ptr<std::vector<double>>& responsibilities);

public:
    constexpr static double CONVERGENCE_THRESHOLD = 0.001;

    // same as above using the default flat prior
    static double logEvidence( const std::vector<std::vector<double>>&  logLikelihoods, double minAF, int nonRefIndex);

    /**
    * @param logLikelihoods matrix of alleles x reads (NOTE: NON_REF allele is assumed to be last)
    * @param priorPseudocounts
    * @param alleleFractionThreshold lower bound of allele fractions to consider for non-ref likelihood
    */
    static double logEvidence(const std::vector<std::vector<double>>&  logLikelihoods, std::vector<double>& priorPseudocounts, double alleleFractionThreshold,  int nonRefIndex);

    /**
    * Given a likelihoods matrix, calculate the parameters of the Dirichlet posterior distribution on their allele
    * fractions, which define a discrete distribution.
    * @param logLikelihoods matrix of alleles x reads
    * @param priorPseudocounts
    */
    static std::shared_ptr<std::vector<double>> alleleFractionsPosterior(const std::vector<std::vector<double>>&  logLikelihoods, std::vector<double>& priorPseudocounts);

    /**
     * Given data log likelihoods and a Dirichlet prior for a categorical distribution, obtain the array of total
     * responsibilities for each category
     * @param logLikelihoods
     * @param dirichletPrior
     * @return
     */
    static std::shared_ptr<std::vector<double>> getEffectiveCounts(const std::vector<std::vector<double>>&  logLikelihoods, std::shared_ptr<std::vector<double>> dirichletPrior);

    static double logDirichletNormalization(std::vector<double>& dirichletParams);

    static std::shared_ptr<std::vector<double>> getColumnOfLogLikelihood(const std::vector<std::vector<double>> &logLikelihoods, int colomnIndex);
};


#endif //MUTECT2CPP_MASTER_SOMATICLIKELIHOODSENGINE_H
