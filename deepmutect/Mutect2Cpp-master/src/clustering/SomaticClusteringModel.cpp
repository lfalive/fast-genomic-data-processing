//
// Created by cluster on 22-11-8.
//

#include "SomaticClusteringModel.h"
#include "NaturalLogUtils.h"
#include "boost/math/distributions.hpp"
#include "BinomialCluster.h"



BetaBinomialCluster SomaticClusteringModel::NEW_CLUSTER(BetaDistributionShape::FLAT_BETA);
BetaDistributionShape SomaticClusteringModel::INITIAL_HIGH_AF_BETA = BetaDistributionShape(10, 1);
BetaDistributionShape SomaticClusteringModel::INITIAL_BACKGROUND_BETA = BetaDistributionShape::FLAT_BETA;
thread_local boost::random::mt19937 SomaticClusteringModel::seed(47382911);

std::vector<double> SomaticClusteringModel::clusterProbabilities(Datum datum) {
    double logVariantPrior = getLogPriorOfSomaticVariant(datum.getIndelLength());
    double logNoVariantPrior = NaturalLogUtils::log1mexp(logVariantPrior);

    std::vector<double> logClusterPosteriors = std::vector<double>(clusters.size()+1, 0);
    for(int i = 0; i < clusters.size() + 1; i++) {
        double logLikelihood = i < clusters.size() ? clusters[i]->logLikelihood(datum) : NEW_CLUSTER.logLikelihood(datum);
        if(i == SEQUENCING_ERROR_INDEX) {
            logClusterPosteriors[i] = logNoVariantPrior + logLikelihood;
        } else if (i == HIGH_AF_INDEX) {
            logClusterPosteriors[i] =  logVariantPrior + logHighAFWeight + logLikelihood;
        } else if (i == BACKGROUND_INDEX) {
            logClusterPosteriors[i] = logVariantPrior + logBackgroundWeight + logLikelihood;
        } else if (i < clusters.size()) {   // existing sparse cluster
            logClusterPosteriors[i] = logVariantPrior + logSparseClustersWeight + logCRPWeight(i)
                   + logLikelihood;
        } else {    // new sparse cluster
            logClusterPosteriors[i] = logVariantPrior + logSparseClustersWeight + logCRPWeight(i)
                   + logLikelihood;
        }
    }

    return  NaturalLogUtils::normalizeLog(logClusterPosteriors, false, false);
}

double SomaticClusteringModel::getLogPriorOfSomaticVariant(int indelLength) {
    if(logVariantPriors.find(indelLength) == logVariantPriors.end()) {
        double input = (*logVariantPriors.begin()).second;
        for(auto & kv : logVariantPriors) {
            if(kv.second < input) {
                input = kv.second;
            }
        }
        logVariantPriors.insert({indelLength, input});
    }
    return logVariantPriors[indelLength] + (indelLength == 0 ? std::log(1.0/3) : 0);
}

SomaticClusteringModel::SomaticClusteringModel(M2FiltersArgumentCollection &MTFAC) : rng(0.0, 1.0){
    logVariantVsArtifactPrior = MTFAC.initialLogPriorOfVariantVersusArtifact;
    totalSparseClusterCount = 0;
    REGULARIZING_PSEUDOCOUNT = 1;
    firstPass = true;
    logHighAFWeight = std::log(INITIAL_HIGH_AF_WEIGHT);
    logBackgroundWeight = std::log(INITIAL_BACKGROUND_WEIGHT);
    std::vector<double> inputs = {logHighAFWeight, logBackgroundWeight};
    logSparseClustersWeight = NaturalLogUtils::log1mexp(NaturalLogUtils::logSumExp(inputs));
    for(int i = -MAX_INDEL_SIZE_IN_PRIOR_MAP; i < MAX_INDEL_SIZE_IN_PRIOR_MAP+1; i++) {
        logVariantPriors.insert({i, MTFAC.getLogIndelPrior()});
    }
    logVariantPriors[0] = MTFAC.getLogSnvPrior();
    clusters = std::vector<std::shared_ptr<AlleleFractionCluster>>(3);
    clusters[SEQUENCING_ERROR_INDEX] = std::shared_ptr<AlleleFractionCluster>(new SequencingError()) ;
    clusters[HIGH_AF_INDEX] = std::shared_ptr<AlleleFractionCluster>(new BetaBinomialCluster(INITIAL_HIGH_AF_BETA));
    clusters[BACKGROUND_INDEX] = std::shared_ptr<AlleleFractionCluster>(new BetaBinomialCluster(INITIAL_BACKGROUND_BETA));
}

double SomaticClusteringModel::logCRPWeight(int clusterIndex) {
    if(clusterIndex < OFFSET) {
        throw std::invalid_argument("Chinese restaurant process does not apply to error, high-AF, and backgorund clusters");
    }
    double numerator = clusterIndex == clusters.size() ? CONCENTRATION : clusterCounts[clusterIndex];
    double denominator = totalSparseClusterCount + CONCENTRATION;
    return std::log(numerator/ denominator);
}

double SomaticClusteringModel::probabilityOfSequencingError(const Datum &datum) {
    return clusterProbabilities(datum)[SEQUENCING_ERROR_INDEX];
}

int SomaticClusteringModel::indelLength(const std::shared_ptr<VariantContext> &vc, int altIndex) {
    return vc->getAlternateAllele(altIndex)->getLength() - vc->getReference()->getLength();
}

SomaticClusteringModel::~SomaticClusteringModel() {

}

void SomaticClusteringModel::record(const std::vector<int> &tumorADs, const std::vector<double> &tumorLogOdds,
                                    double artifactProbability, double nonSomaticProbability,
                                    const std::shared_ptr<VariantContext> &vc) {
    int totalAD = 0;
    for(int i = 0; i < tumorADs.size(); i++) {
        totalAD += tumorADs[i];
    }
    for(int i = 0; i < tumorLogOdds.size(); i++) {
        data.emplace_back(tumorLogOdds[i], artifactProbability, nonSomaticProbability, tumorADs[i+1], totalAD,
                           indelLength(vc, i));
    }
}

void SomaticClusteringModel::learnAndClearAccumulatedData() {
    if(firstPass) {
        clusterAssignments = std::vector<std::optional<int>>(data.size(), std::nullopt);
        for(int i = 0; i < clusters.size(); i++) {
            clusterCounts.emplace_back(0);
        }
    }
    for(int iteration = 0; iteration < NUM_ITERATIONS; iteration++) {
        for (int datumIndex = 0; datumIndex < data.size(); datumIndex++) {
            Datum datum = popDatum(datumIndex);
            if(rng(seed) < datum.getNonSequencingErrorProb()) {
                continue;
            }
            std::vector<double> clusterPosteriors = clusterProbabilities(datum);
            std::uniform_real_distribution<double> urd(0, 1);
            double tmp = urd(seed);
            int clusterIndex = 0;
            double sum = 0;
            while(sum < tmp) {
                sum += clusterPosteriors[clusterIndex];
                clusterIndex++;
            }
            assignDatum(datumIndex, --clusterIndex);
        }
        pruneEmptyClusters();
        std::vector<std::vector<Datum>> dataByCluster = std::vector<std::vector<Datum>>(clusters.size());
        for(int i = 0; i < clusterAssignments.size(); i++) {
            if(clusterAssignments[i].has_value()) {
                dataByCluster[clusterAssignments[i].value()].emplace_back(data[i]);
            }
        }
        for(int i = 0; i < clusters.size(); i++) {
            clusters[i]->learn(dataByCluster[i]);
        }
        learnWeightsAndPriors();
    }

    firstPass = false;
    data.clear();
}

Datum SomaticClusteringModel::popDatum(int datumIndex) {
    if(clusterAssignments[datumIndex].has_value()) {
        int c = clusterAssignments[datumIndex].value();
        clusterCounts[c]--;
        if(OFFSET <= c) {
            totalSparseClusterCount--;
        }
    }
    clusterAssignments[datumIndex] = std::nullopt;
    return data[datumIndex];
}

void SomaticClusteringModel::assignDatum(int datumIndex, int clusterIndex) {
    Datum datum = data[datumIndex];
    if(clusterIndex == clusters.size()) {
        float rand_uniform = rng(seed);

        boost::math::beta_distribution<> beta_dist(datum.getAltCount()+1, datum.getTotalCount()-datum.getAltCount()+1);
        auto newClusterAlleleFraction = boost::math::quantile(beta_dist, rand_uniform);
        clusters.emplace_back(new BinomialCluster(newClusterAlleleFraction));
        clusterCounts.emplace_back(0);
    }

    if(OFFSET <= clusterIndex) {
        totalSparseClusterCount++;
    }

    clusterAssignments[datumIndex] = clusterIndex;
    clusterCounts[clusterIndex]++;
}

void SomaticClusteringModel::pruneEmptyClusters() {
    std::map<int, int> oldToNewClusterIndices;
    for(int i = 0; i < OFFSET; i++) {
        oldToNewClusterIndices.insert({i, i});
    }
    int newIndex = OFFSET;
    for (int oldIndex = OFFSET; oldIndex < clusters.size(); oldIndex++) {
        if (clusterCounts[oldIndex] > 0) {
            oldToNewClusterIndices[oldIndex] = newIndex;

            if (newIndex != oldIndex) {
                clusters[newIndex] = clusters[oldIndex];
                clusterCounts[newIndex] = clusterCounts[oldIndex];
            }
            newIndex++;
        }
    }

    int iter = std::min(newIndex, (int)clusters.size());
    std::vector<std::shared_ptr<AlleleFractionCluster>> tmp1;
    for(int i = 0; i < newIndex; i++) {
        tmp1.emplace_back(clusters[i]);
    }
    std::swap(tmp1, clusters);
    std::vector<int> tmp2;
    for(int i = 0; i < newIndex; i++) {
        tmp2.emplace_back(clusterCounts[i]);
    }
    std::swap(tmp2, clusterCounts);
    std::vector<std::optional<int>> tmp;
    for(auto & k : clusterAssignments) {
        if(k.has_value()) {
            tmp.emplace_back(oldToNewClusterIndices[k.value()]);
        } else {
            tmp.emplace_back(k);
        }
    }
    std::swap(clusterAssignments, tmp);
}

void SomaticClusteringModel::learnWeightsAndPriors() {
    double totalVariants = clusterCounts[HIGH_AF_INDEX] + clusterCounts[BACKGROUND_INDEX] + totalSparseClusterCount
            + REGULARIZING_PSEUDOCOUNT;
    logHighAFWeight = std::log(REGULARIZING_PSEUDOCOUNT + clusterCounts[HIGH_AF_INDEX] / totalVariants);
    logBackgroundWeight = std::log((REGULARIZING_PSEUDOCOUNT + clusterCounts[BACKGROUND_INDEX]) / totalVariants);
    logSparseClustersWeight = std::log((REGULARIZING_PSEUDOCOUNT + totalSparseClusterCount) /totalVariants);

    std::vector<int> tmpLength;
    for(int i = 0; i < data.size(); i++) {
        if(clusterAssignments[i].value_or(0) != 0) {
            tmpLength.emplace_back(data[i].getIndelLength());
        }
    }
    std::map<int, long> variantCountsByIndelLength;
    for(int length : tmpLength) {
        variantCountsByIndelLength[length]++;
    }
    double technicalArtifactCount = 0;
    for(auto & da : data) {
        technicalArtifactCount += da.getArtifactProb();
    }
    if(callableSites.has_value()) {
        for(int i = -MAX_INDEL_SIZE_IN_PRIOR_MAP; i <= MAX_INDEL_SIZE_IN_PRIOR_MAP; i++) {
            long value = variantCountsByIndelLength.find(i) != variantCountsByIndelLength.end() ? variantCountsByIndelLength[i] : 0;
            double empiricalRatio = value / callableSites.value();
            logVariantPriors[i] = std::log(std::max(empiricalRatio, i == 0 ? 1.0e-8 : 1.0e-9));
        }
    }
    long variantCount = 0;
    for(auto & k : variantCountsByIndelLength) {
        variantCount += k.second;
    }
    logVariantVsArtifactPrior = std::log((variantCount + REGULARIZING_PSEUDOCOUNT) / (variantCount + technicalArtifactCount + REGULARIZING_PSEUDOCOUNT * 2));
}

double SomaticClusteringModel::logLikelihoodGivenSomatic(int totalCount, int altCount) {
    std::vector<double> logClusterLikelihoods;
    for(int i = 0; i < clusters.size(); i++) {
        if(i != SEQUENCING_ERROR_INDEX) {
            double logLikelihood = clusters[i]->logLikelihood(totalCount, altCount);
            if (i == HIGH_AF_INDEX) {
                logClusterLikelihoods.emplace_back(logHighAFWeight + logLikelihood);
            } else if (i == BACKGROUND_INDEX) {
                logClusterLikelihoods.emplace_back(logBackgroundWeight + logLikelihood);
            } else {   // sparse cluster
                logClusterLikelihoods.emplace_back(logSparseClustersWeight + logCRPWeight(i) + logLikelihood);
            }
        }
    }
    return NaturalLogUtils::logSumExp(logClusterLikelihoods);
}

double SomaticClusteringModel::getLogPriorOfSomaticVariant(const std::shared_ptr<VariantContext> &vc, int altIndex) {
    int indellen = indelLength(vc, altIndex);
    return getLogPriorOfSomaticVariant(indellen);
}
