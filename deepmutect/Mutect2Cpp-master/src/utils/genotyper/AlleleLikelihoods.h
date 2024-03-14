//
// Created by lhh on 4/23/22.
//

#ifndef MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H
#define MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H

#include <utility>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <cassert>
#include <functional>
#include "parallel_hashmap/phmap.h"
#include <cassert>
#include "Allele.h"
#include "samtools/SAMRecord.h"
#include "MathUtils.h"
#include "Fragment.h"
#include "Haplotype.h"
#include "parallel_hashmap/phmap.h"

using namespace std;

template <typename E, typename A>
class SampleMatrix;
template <typename E, typename A>
class BestAllele;

template <typename E, typename A>
class AlleleLikelihoods {
    friend class SampleMatrix<E, A>;
    friend class BestAllele<E, A>;

private:
    const static int MISSING_REF = -1;

    /**
    * Index of the reference allele if any, otherwise {@link #MISSING_REF}.
    */
    int referenceAlleleIndex = MISSING_REF;

    /**
     * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
     */
    vector<SampleMatrix<E, A>*> sampleMatrices;



    double getInformativeThreshold() {
        return isNaturalLog ? NATURAL_LOG_INFORMATIVE_THRESHOLD : LOG_10_INFORMATIVE_THRESHOLD;
    }

    // Search for the reference allele, if not found the index is {@link MISSING_REF}.
    static int findReferenceAllele(std::vector<shared_ptr<A>>& alleles){
        int number = alleles.size();
        for(int i=0; i<number; i++)
        {
            if(alleles[i]->getIsReference())
                return i;
        }
        return MISSING_REF;
    }

    void setupIndexes(std::map<std::string, std::vector<std::shared_ptr<E>>>& evidenceBySample, int sampleCount, int alleleCount) {
        for(int s = 0; s < sampleCount; s++)
        {
            std::string& sample = samples[s];

            evidenceBySampleIndex->emplace_back(evidenceBySample[sample]);

            int sampleEvidenceCount = (*evidenceBySampleIndex)[s].size();

            valuesBySampleIndex->emplace_back(vector<vector<double>>(alleleCount, vector<double>(sampleEvidenceCount, 0.0)));
        }
    }



    /**
     * Search the best allele for a unit of evidence.
     *
     * @param sampleIndex including sample index.
     * @param evidenceIndex  target evidence index.
     *
     * @param priorities An array of allele priorities (higher values have higher priority) to be used, if present, to break ties for
     *                   uninformative likelihoods, in which case the evidence is assigned to the allele with the higher score.
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null}
     * if non-could be found.
     */
    shared_ptr<BestAllele<E, A>> searchBestAllele(int sampleIndex, int evidenceIndex, bool canBeReference, const double* priorities){
        int alleleCount = alleles->size();
        if (alleleCount == 0 || (alleleCount == 1 && referenceAlleleIndex == 0 && !canBeReference)) {
            return make_shared<BestAllele<E, A>>(this, sampleIndex, evidenceIndex, -1, -numeric_limits<double>::infinity(), -numeric_limits<double>::infinity());
        }

        auto & sampleValues = (*valuesBySampleIndex)[sampleIndex];
        int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;

        int secondBestIndex = 0;
        double bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        double secondBestLikelihood =  -numeric_limits<double>::infinity();

        for(int a = bestAlleleIndex + 1; a < alleleCount; a++){
            if (!canBeReference && referenceAlleleIndex == a) {
                continue;
            }
            double candidateLikelihood = sampleValues[a][evidenceIndex];
            if (candidateLikelihood > bestLikelihood) {
                secondBestIndex = bestAlleleIndex;
                bestAlleleIndex = a;
                secondBestLikelihood = bestLikelihood;
                bestLikelihood = candidateLikelihood;
            } else if (candidateLikelihood > secondBestLikelihood) {
                secondBestIndex = a;
                secondBestLikelihood = candidateLikelihood;
            }
        }

        if (priorities != nullptr && bestLikelihood - secondBestLikelihood < getInformativeThreshold()) {
            double bestPriority = priorities[bestAlleleIndex];
            double secondBestPriority = priorities[secondBestIndex];
            for (int a = 0; a < alleleCount; a++) {
                double candidateLikelihood = sampleValues[a][evidenceIndex];
                if (a == bestAlleleIndex || (!canBeReference && a == referenceAlleleIndex) || bestLikelihood - candidateLikelihood > getInformativeThreshold()) {
                    continue;
                }
                double candidatePriority = priorities[a];

                if (candidatePriority > bestPriority) {
                    secondBestIndex = bestAlleleIndex;
                    bestAlleleIndex = a;
                    secondBestPriority = bestPriority;
                    bestPriority = candidatePriority;
                } else if (candidatePriority > secondBestPriority) {
                    secondBestIndex = a;
                    secondBestPriority = candidatePriority;
                }
            }
        }

        bestLikelihood = sampleValues[bestAlleleIndex][evidenceIndex];
        secondBestLikelihood = secondBestIndex != bestAlleleIndex ? sampleValues[secondBestIndex][evidenceIndex] : -numeric_limits<double>::infinity();
        return make_shared<BestAllele<E, A>>(this, sampleIndex, evidenceIndex, bestAlleleIndex, bestLikelihood, secondBestLikelihood);
    }

    shared_ptr<BestAllele<E, A>> searchBestAllele(int sampleIndex, int evidenceIndex, bool canBeReference){
        return searchBestAllele(sampleIndex, evidenceIndex, canBeReference, nullptr);
    }

    // Does the normalizeLikelihoods job for each piece of evidence.
    void normalizeLikelihoodsPerEvidence(double maximumBestAltLikelihoodDifference, vector<vector<double>>& sampleValues, int sampleIndex, int evidenceIndex){
        //allow the best allele to be the reference because asymmetry leads to strange artifacts like het calls with >90% alt reads
        shared_ptr<BestAllele<E, A>> bestAllele = searchBestAllele(sampleIndex, evidenceIndex, true);

        double worstLikelihoodCap = bestAllele->likelihood + maximumBestAltLikelihoodDifference;

        int alleleCount = alleles->size();

        // Guarantee to be the case by enclosing code.
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] < worstLikelihoodCap) {
                sampleValues[a][evidenceIndex] = worstLikelihoodCap;
            }
        }
    }

    phmap::flat_hash_map<shared_ptr<E>, int>& getEvidenceIndexBySampleIndex(int sampleIndex)
    {
        if(evidenceIndexBySampleIndex.size() <= sampleIndex)
        {
            for(int i = evidenceIndexBySampleIndex.size()-1; i<sampleIndex+1; i++)
                evidenceIndexBySampleIndex.push_back({});
        }

        if(evidenceIndexBySampleIndex[sampleIndex].empty()){
            auto& sampleEvidence = evidenceBySampleIndex->operator[](sampleIndex);
            int sampleEvidenceCount = sampleEvidence.size();
            for(int r=0; r<sampleEvidenceCount; r++)
            {
                evidenceIndexBySampleIndex[sampleIndex].template emplace(sampleEvidence[r], r);
            }
        }
        return evidenceIndexBySampleIndex[sampleIndex];
    }

protected:
    bool isNaturalLog = false;

    /**
     * Evidence by sample index. Each sub array contains reference to the evidence of the ith sample.
     */
    shared_ptr<std::vector<std::vector<std::shared_ptr<E>>>> evidenceBySampleIndex;

    shared_ptr<std::vector<std::vector<std::vector<double>>>> valuesBySampleIndex;

    std::vector<std::string> samples;

    shared_ptr<vector<shared_ptr<A>>> alleles;

    /**
    * Maps from each unit of evidence to its index within the sample.
    *
    * <p>In order to save CPU time the indices contained in this array (not the array itself) is
    * lazily initialized by invoking {@link #evidenceIndexBySampleIndex(int)}.</p>
    */
    vector<phmap::flat_hash_map<shared_ptr<E>, int>> evidenceIndexBySampleIndex;

    double maximumLikelihoodOverAllAlleles(int sampleIndex, int evidenceIndex) {
        double result = -numeric_limits<double>::infinity();
        int alleleCount = alleles->size();
        auto & sampleValues = (*valuesBySampleIndex)[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][evidenceIndex] > result) {
                result = sampleValues[a][evidenceIndex];
            }
        }
        return result;
    }

public:
    constexpr static double LOG_10_INFORMATIVE_THRESHOLD = 0.2;
    static double NATURAL_LOG_INFORMATIVE_THRESHOLD;

    AlleleLikelihoods(vector<string>& samples, shared_ptr<vector<shared_ptr<A>>>& alleles, map<string, vector<std::shared_ptr<E>>>& evidenceBySample) :
    samples(samples), alleles(alleles), evidenceBySampleIndex(new std::vector<std::vector<std::shared_ptr<E>>>())
    {
        int sampleCount = samples.size();
        int alleleCount = alleles->size();

        evidenceBySampleIndex->reserve(sampleCount);
        valuesBySampleIndex = make_shared<vector<vector<vector<double>>>>();
        valuesBySampleIndex->reserve(sampleCount);
        referenceAlleleIndex = findReferenceAllele(*alleles);

        setupIndexes(evidenceBySample, sampleCount, alleleCount);
        sampleMatrices.reserve(sampleCount);
        for(int i=0; i<sampleCount; i++)
        {
            sampleMatrices.template emplace_back(nullptr);
        }

    }

    AlleleLikelihoods(std::shared_ptr<std::vector<shared_ptr<A>>> alleles, std::vector<std::string>& samples, shared_ptr<vector<vector<shared_ptr<E>>>> evidenceBySampleIndex, shared_ptr<vector<vector<vector<double>>>> values)
        : samples(samples), alleles(alleles), evidenceBySampleIndex(evidenceBySampleIndex), valuesBySampleIndex(std::move(values))
    {
        int sampleCount = samples.size();
        evidenceIndexBySampleIndex.reserve(sampleCount);

        referenceAlleleIndex = findReferenceAllele(*alleles);
        sampleMatrices.reserve(sampleCount);
        for(int i=0; i<sampleCount; i++)
        {
            sampleMatrices.template emplace_back(nullptr);
        }
    }

    ~AlleleLikelihoods()
    {
        for(auto sampleMatrix : sampleMatrices)
        {
            delete sampleMatrix;
        }
    }

    /**
     * Returns the units of evidence that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex the requested sample.
     * @return never {@code null} but perhaps a zero-length array if there is no evidence in sample. No element in
     *   the array will be null.
     */
    vector<shared_ptr<E>>& sampleEvidence(int sampleIndex){
        return (*evidenceBySampleIndex)[sampleIndex];
    }


    /**
     * Returns an evidence vs allele likelihood matrix corresponding to a sample.
     */
    SampleMatrix<E, A> * sampleMatrix(int sampleIndex)
    {
        assert(sampleIndex >= 0 && sampleIndex < samples.size());

        if(!sampleMatrices[sampleIndex])
            sampleMatrices[sampleIndex] = new SampleMatrix<E, A>(sampleIndex, this);
        return sampleMatrices[sampleIndex];
    }

    /**
     * Returns the samples in this evidence-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable. It will not be updated if the collection
     *     allele list changes.
     * </p>
     *
     * @return never {@code null}.
     */
    vector<shared_ptr<A>>& getAlleles(){
        return *alleles;
    }

    shared_ptr<vector<shared_ptr<A>>> getAlleleList()
    {
        return alleles;
    }

    /**
    * Returns the allele given its index.
    *
    * @param alleleIndex the allele index.
    *
    * @throws IllegalArgumentException the allele index is {@code null}.
    *
    * @return never {@code null}.
    */
    shared_ptr<A> getAllele(int alleleIndex)
    {
        assert(alleleIndex >= 0 && alleleIndex < alleles->size());
        return alleles->operator[](alleleIndex);
    }

    void setIsNaturalLog(bool isNaturalLog)
    {
        this->isNaturalLog = isNaturalLog;
    }

    int indexOfSamples(string& sample)
    {
        for(int i=0; i<samples.size(); i++)
        {
            if(sample == samples[i])
                return i;
        }
        throw "Sample not found";
    }

    /**
    * Number of samples included in the likelihood collection.
    * @return 0 or greater.
    */
    int numberOfSamples()
    {
        return samples.size();
    }

    int numberOfAlleles()
    {
        return alleles->size();
    }

    string & getSample(int sampleIndex)
    {
        return samples[sampleIndex];
    }

    /**
     * Adjusts likelihoods so that for each unit of evidence, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each unit of evidence based on the maximum alternative allele likelihood.
     *
     * @param maximumLikelihoodDifferenceCap maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    void normalizeLikelihoods(double maximumLikelihoodDifferenceCap){
        assert(!isnan(maximumLikelihoodDifferenceCap) && maximumLikelihoodDifferenceCap < 0.0);

        if(isinf(maximumLikelihoodDifferenceCap))
            return;

        int alleleCount = alleles->size();
        if(alleleCount == 0 || alleleCount == 1)
            return;

        for(int s=0; s<valuesBySampleIndex->size(); s++){
            auto& sampleValues = (*valuesBySampleIndex)[s];
            int evidenceCount = (*evidenceBySampleIndex)[s].size();
            for(int r=0; r<evidenceCount; r++)
            {
                normalizeLikelihoodsPerEvidence(maximumLikelihoodDifferenceCap, sampleValues, s, r);
            }
        }
    }

    /**
   * Removes those read that the best possible likelihood given any allele is just too low.
   *
   * <p>
   *     This is determined by a maximum error per read-base against the best likelihood possible.
   * </p>
   *
   * @param log10MinTrueLikelihood Function that returns the minimum likelihood that the best allele for a unit of evidence must have
   * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
   *
   * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
   */   // TODO: validate this method
    void filterPoorlyModeledEvidence(const function<double(shared_ptr<SAMRecord>, double)>& log10MinTrueLikelihood, double maximumErrorPerBase){
        assert(alleles->size() > 0);
        int numberOfSamples = samples.size();
        for (int s = 0; s < numberOfSamples; s++) {
            auto & sampleEvidence = (*evidenceBySampleIndex)[s];
            vector<int> indexesToRemove;

            int numberOfEvidence = sampleEvidence.size();
            for (int e = 0; e < numberOfEvidence; e++) {
                if (maximumLikelihoodOverAllAlleles(s, e) <  log10MinTrueLikelihood(sampleEvidence[e], maximumErrorPerBase)){
                    indexesToRemove.push_back(e);
                }
            }
            removeSampleEvidence(s, indexesToRemove, alleles->size());
        }
    }

    // ---we use indexes to be removed instead of evidences to be removed
    // Requires that the collection passed iterator can remove elements, and it can be modified.
    void removeSampleEvidence(int sampleIndex, vector<int>& indexesToRemove, int alleleCount) {
        if(indexesToRemove.empty())
            return;

        auto& sampleEvidence = (*evidenceBySampleIndex)[sampleIndex];

        vector<int> IndicesToKeep;
        int numberOfEvidence = sampleEvidence.size();
        int index = 0;
        for(int i=0; i<numberOfEvidence; i++)
        {
            if(index >= indexesToRemove.size() || indexesToRemove[index] != i)
            {
                IndicesToKeep.push_back(i);
            } else {
                index++;
            }
        }

        // Then we skim out the likelihoods of the removed evidence.
        auto & oldSampleValues = (*valuesBySampleIndex)[sampleIndex];
        vector<vector<double>> newSampleValues(alleleCount, vector<double>(IndicesToKeep.size(), 0.0));

        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < IndicesToKeep.size(); r++) {
                newSampleValues[a][r] = oldSampleValues[a][IndicesToKeep[r]];
            }
        }

        vector<shared_ptr<E>> newSampleEvidence;
        for(int & i : IndicesToKeep)
        {
            newSampleEvidence.template emplace_back(sampleEvidence[i]);
        }

        (*valuesBySampleIndex)[sampleIndex] = newSampleValues;
        (*evidenceBySampleIndex)[sampleIndex] = newSampleEvidence;
    }

    void switchToNaturalLog(){
        assert(!isNaturalLog);
        int sampleCount = samples.size();
        int alleleCount = alleles->size();

        for (int s = 0; s < sampleCount; s++) {
            int evidenceCount = sampleEvidenceCount(s);
            for (int a = 0; a < alleleCount; a++) {
                for (int e = 0; e < evidenceCount; e++) {
                    (*valuesBySampleIndex)[s][a][e] = MathUtils::log10ToLog((*valuesBySampleIndex)[s][a][e]);
                }
            }
        }
        isNaturalLog = true;
    }

    /**
    * Returns the quantity of evidence that belongs to a sample in the evidence-likelihood collection.
    * @param sampleIndex the query sample index.
    *
    * @return 0 or greater.
    */
    int sampleEvidenceCount(int sampleIndex) {
        assert(sampleIndex >= 0 && sampleIndex < samples.size());
        return (*evidenceBySampleIndex)[sampleIndex].size();
    }

    /**
     * Default version where ties are broken in favor of the reference allele
     */
    shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> bestAllelesBreakingTies()
    {
        auto result = make_shared<vector<shared_ptr<BestAllele<E, A>>>>();
        for(int i=0; i<numberOfSamples(); i++)
        {
            auto bestAlleles = bestAllelesBreakingTies(i);
            for(auto& bestAllele: *bestAlleles)
            {
                result->template emplace_back(bestAllele);
            }
        }
        return result;
    }

    shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> bestAllelesBreakingTies(std::string sample)
    {
        int sampleIndex= indexOfSamples(sample);
        return bestAllelesBreakingTies(sampleIndex);
    }

    shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> bestAllelesBreakingTies(int sampleIndex)
    {
        return bestAllelesBreakingTies(sampleIndex, [](shared_ptr<A> a){return a->getIsReference() ? 1.0 : 0.0;});
    }

    /**
    * Returns the collection of best allele estimates for the evidence based on the evidence-likelihoods.
    * "Ties" where the ref likelihood is within {@code AlleleLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
    * are broken by the {@code tieBreakingPriority} function.
    *
    * @return never {@code null}, one element per unit fo evidence in the evidence-likelihoods collection.
    */
    shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> bestAllelesBreakingTies(function<double(shared_ptr<A>)> tieBreakingPriority)
    {
        vector<double> priorities;
        if(!alleles->empty())
        {
            for(shared_ptr<A> h : *alleles)
            {
                assert(h != nullptr);
                double tmp = tieBreakingPriority(h);
                priorities.push_back(tmp);

            }

        }

        shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> result = shared_ptr<vector<shared_ptr<BestAllele<E, A>>>>(new vector<shared_ptr<BestAllele<E, A>>>);
        for(int sampleIndex = 0; sampleIndex < samples.size(); sampleIndex++)
        {
            int evidenceCount = (*evidenceBySampleIndex)[sampleIndex].size();

            for(int r=0; r<evidenceCount; r++)
            {
                result->template emplace_back(searchBestAllele(sampleIndex, r, true, priorities.data()));
            }
        }

        return result;
    }

    shared_ptr<vector<shared_ptr<BestAllele<E, A>>>> bestAllelesBreakingTies(int sampleIndex, function<double(shared_ptr<A>)> tieBreakingPriority)
    {
        //TODO: this currently just does ref vs alt.  Really we want CIGAR complexity.
        double * priorities = nullptr;
        if(alleles != nullptr)
        {
            priorities = new double[alleles->size()];
            for(int i=0; i<alleles->size(); i++)
            {
                priorities[i] = tieBreakingPriority(alleles->operator[](i));
            }
        }

        int evidenceCount = evidenceBySampleIndex->operator[](sampleIndex).size();
        auto result = make_shared<vector<shared_ptr<BestAllele<E, A>>>>();
        result->reserve(evidenceCount);
        for(int r=0; r<evidenceCount; r++)
            result->template emplace_back(searchBestAllele(sampleIndex, r, true, priorities));

        delete[] priorities;
        return result;
    }

    // TODO: is evidenceIndexBySampleIndex useful ?
    void changeEvidence(shared_ptr<phmap::flat_hash_map<shared_ptr<E>, shared_ptr<E>>> evidenceReplacements)
    {
        int sampleCount = samples.size();
        for(int s = 0; s < sampleCount; s++)
        {
            auto & sampleEvidence = (*evidenceBySampleIndex)[s];
            // Object2IntMap<EVIDENCE> evidenceIndex = evidenceIndexBySampleIndex.get(s);
            int sampleEvidenceCount = sampleEvidence.size();
            for (int r = 0; r < sampleEvidenceCount; r++) {
                shared_ptr<E>& evidence = sampleEvidence[r];
                shared_ptr<E>& replacement = evidenceReplacements->at(evidence);
                if(replacement == nullptr)
                    continue;

                sampleEvidence[r] = replacement;
               /* if (evidenceIndex != null) {
                    evidenceIndex.remove(evidence);
                    evidenceIndex.put(replacement, r);
                }*/
            }
        }
    }

    int evidenceIndex(int sampleIndex, shared_ptr<E>& evidence)
    {
        /*if(evidenceIndexBySampleIndex[sampleIndex].template find(evidence) != evidenceIndexBySampleIndex[sampleIndex].end())
        {
            return evidenceIndexBySampleIndex[sampleIndex].at(evidence);
        }
        return -1;*/
        auto& map = getEvidenceIndexBySampleIndex(sampleIndex);
        if(map.template find(evidence) != map.end())
            return map.at(evidence);
        return -1;
    }

    /**
     * Group evidence into lists of evidence -- for example group by read name to force read pairs to support a single haplotype.
     *
     * Log Likelihoods are summed over all evidence
     * in a group, corresponding to an independent evidence assumption.  Since this container's likelihoods generally pertain to
     * sequencing only (and not sample prep etc) this is usually a good assumption.
     *
     * @param groupingFunction Attribute function for grouping evidence, for example GATKRead::getName
     * @param gather Transformation applied to collections of evidence with same value of groupingFunction.  For example, Fragment::new
     *               to construct a fragment out of a pair of reads with the same name
     *
     * @return a new AlleleLikelihoods based on the grouped, transformed evidence.
    */
    AlleleLikelihoods<Fragment, A>* groupEvidence(function<std::string&(shared_ptr<E>)> groupingFunction, function<shared_ptr<Fragment>(vector<shared_ptr<E>>&)> gather)
     {
         int sampleCount = samples.size();
         auto newLikelihoodValues = make_shared<vector<vector<vector<double>>>>();
         int alleleCount = alleles->size();

         auto newEvidenceBySampleIndex = make_shared<vector<vector<shared_ptr<Fragment>>>>(sampleCount, vector<shared_ptr<Fragment>>());

         for(int s = 0; s < sampleCount; s++)
         {
             vector<shared_ptr<vector<shared_ptr<E>>>> evidenceGroups;  //---Maybe this variable is unnecessary
             vector<shared_ptr<E>> & sampleEvidence = (*evidenceBySampleIndex)[s];
             phmap::flat_hash_map<string, shared_ptr<vector<shared_ptr<E>>>> map;
             for(auto& evidence : sampleEvidence)
             {
                 string & groupingKey = groupingFunction(evidence);
                 if(map.find(groupingKey) != map.end())
                 {
                     map[groupingKey]->emplace_back(evidence);
                 } else {
                     map.template emplace(groupingKey, make_shared<vector<shared_ptr<E>>>(1, evidence));
                 }

             }

             for(auto& kv: map)
                 evidenceGroups.push_back(std::move(kv.second));


             int newEvidenceCount = evidenceGroups.size();
             auto& oldSampleValues = (*valuesBySampleIndex)[s];
             newLikelihoodValues->template emplace_back(vector<vector<double>>(alleleCount, vector<double>(newEvidenceCount, 0.0)));

             // For each old allele and read we update the new table keeping the maximum likelihood.
             for (int a = 0; a < alleleCount; a++) {
                 for (int newEvidenceIndex = 0; newEvidenceIndex < newEvidenceCount; newEvidenceIndex++) {
                     for(auto& evidence : *evidenceGroups[newEvidenceIndex])
                     {
                         int oldEvidenceIndex = evidenceIndex(s, evidence);
                         assert(oldEvidenceIndex != -1);
                         newLikelihoodValues->operator[](s)[a][newEvidenceIndex] += oldSampleValues[a][oldEvidenceIndex];
                     }
                 }
             }

             vector<shared_ptr<Fragment>>& temp = newEvidenceBySampleIndex->operator[](s);
             for(auto& group : evidenceGroups)
             {
                 temp.template emplace_back(gather(*group));
             }
         }

         // Finally we create the new read-likelihood
         auto* result = new AlleleLikelihoods<Fragment, A>(alleles, samples, newEvidenceBySampleIndex, newLikelihoodValues);

         return result;
     }

     /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each unit of evidence in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and evidence as the original.
     *
     * @param overlap if not {@code null}, only units of evidence that overlap the location (with unclipping) will be present in
     *                        the output evidence-collection.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this evidence-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
     AlleleLikelihoods<E, Allele>* marginalize(std::shared_ptr<phmap::flat_hash_map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Haplotype>>>, hash_Allele, equal_Allele>> newToOldAlleleMap, shared_ptr<SimpleInterval> overlap)
     {
         assert(newToOldAlleleMap != nullptr);
         if(overlap == nullptr)
             return marginalize(newToOldAlleleMap);

         shared_ptr<vector<shared_ptr<Allele>>> newAlleles(new vector<shared_ptr<Allele>>);
         for(auto & iter : *newToOldAlleleMap)
         {
             newAlleles->template emplace_back(iter.first);
         }
         int oldAlleleCount = alleles->size();
         int newAlleleCount = newAlleles->size();

         // we get the index correspondence between new old -> new allele, -1 entries mean that the old
         // allele does not map to any new; supported but typically not the case.
         auto oldToNewIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, *newAlleles);

         // We calculate the marginal likelihoods.
         auto evidenceToKeep = overlappingEvidenceIndicesBySampleIndex(overlap);
         shared_ptr<vector<vector<vector<double>>>> newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, *oldToNewIndexMap, evidenceToKeep);
         int sampleCount = samples.size();

         shared_ptr<vector<vector<shared_ptr<E>>>> newEvidenceBySampleIndex(new vector<vector<shared_ptr<E>>>());

         for (int s = 0; s < sampleCount; s++) {
             vector<int>& sampleEvidenceToKeep = (*evidenceToKeep)[s];
             vector<shared_ptr<E>>& oldSampleEvidence = (*evidenceBySampleIndex)[s];
             int oldSampleEvidenceCount = oldSampleEvidence.size();
             int newSampleEvidenceCount = sampleEvidenceToKeep.size();

             vector<shared_ptr<E>> newSampleEvidence;
             if(newSampleEvidenceCount == oldSampleEvidenceCount)
             {
                 newEvidenceBySampleIndex->template emplace_back(oldSampleEvidence);
             } else {
                 for(int i=0; i<sampleEvidenceToKeep.size(); i++)
                 {
                     newSampleEvidence.emplace_back(oldSampleEvidence[sampleEvidenceToKeep[i]]);
                 }
                 newEvidenceBySampleIndex->template emplace_back(newSampleEvidence);
             }
         }

         // Finally we create the new evidence-likelihood
         auto * result = new AlleleLikelihoods<E, Allele>(newAlleles, samples, newEvidenceBySampleIndex, newLikelihoodValues);
         result->setIsNaturalLog(isNaturalLog);
         return result;
     }


     AlleleLikelihoods<E, Allele>* marginalize(std::shared_ptr<phmap::flat_hash_map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Haplotype>>>, hash_Allele, equal_Allele>> newToOldAlleleMap)
     {
        assert(newToOldAlleleMap != nullptr);
        shared_ptr<vector<shared_ptr<Allele>>> newAlleles = make_shared<vector<shared_ptr<Allele>>>();
        for(auto & iter : *newToOldAlleleMap)
        {
            newAlleles->emplace_back(iter.first);
        }
        int oldAlleleCount = alleles->size();
        int newAlleleCount = newAlleles->size();

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        auto _oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, *newAlleles);

        // We calculate the marginal likelihoods.
        auto newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, *_oldToNewAlleleIndexMap, nullptr);
        int sampleCount = samples.size();

        auto newEvidenceBySampleIndex = make_shared<vector<vector<shared_ptr<E>>>>(sampleCount, vector<shared_ptr<E>>());

        for (int s = 0; s < sampleCount; s++) {
            newEvidenceBySampleIndex->template emplace_back(evidenceBySampleIndex->operator[](s));
        }

        // Finally we create the new evidence-likelihood
        auto* result = new AlleleLikelihoods<E, Allele>(
                newAlleles,
                samples,
                newEvidenceBySampleIndex,
                newLikelihoodValues);
        result->setIsNaturalLog(isNaturalLog);
        return result;
     }

     // pay attention! The data type of parameters is different from the methods above
     AlleleLikelihoods<E, Allele>* marginalize(shared_ptr<std::map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Allele>>>>> newToOldAlleleMap)
     {
         assert(newToOldAlleleMap != nullptr);
         shared_ptr<vector<shared_ptr<Allele>>> newAlleles = make_shared<vector<shared_ptr<Allele>>>();
         for(auto & iter : *newToOldAlleleMap)
         {
             newAlleles->emplace_back(iter.first);
         }
         int oldAlleleCount = alleles->size();
         int newAlleleCount = newAlleles->size();

         // we get the index correspondence between new old -> new allele, -1 entries mean that the old
         // allele does not map to any new; supported but typically not the case.
         auto _oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, *newAlleles);

         // We calculate the marginal likelihoods.
         auto newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, *_oldToNewAlleleIndexMap, nullptr);
         int sampleCount = samples.size();

         auto newEvidenceBySampleIndex = make_shared<vector<vector<shared_ptr<E>>>>();

         for (int s = 0; s < sampleCount; s++) {
             newEvidenceBySampleIndex->template emplace_back(evidenceBySampleIndex->operator[](s));
         }

         // Finally we create the new evidence-likelihood
         auto* result = new AlleleLikelihoods<E, Allele>(
                 newAlleles,
                 samples,
                 newEvidenceBySampleIndex,
                 newLikelihoodValues);
         result->setIsNaturalLog(isNaturalLog);
         return result;
     }

     // calculates an old to new allele index map array.
     shared_ptr<vector<int>> oldToNewAlleleIndexMap(std::shared_ptr<phmap::flat_hash_map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Haplotype>>>, hash_Allele, equal_Allele>> newToOldAlleleMap, int oldAlleleCount,  vector<shared_ptr<Allele>>& newAlleles)
     {
         for(auto& h : newAlleles)
             assert(h);

         auto oldToNewAlleleIndexMap = make_shared<vector<int>>(oldAlleleCount, -1); // -1 indicate that there is no new allele that make reference to that old one.
         for(int newIndex = 0; newIndex < newAlleles.size(); newIndex++)
         {
            auto newAllele = newAlleles[newIndex];
            for(shared_ptr<Haplotype> oldAllele : *newToOldAlleleMap->at(newAllele))
            {
                int oldAlleleIndex = indexOfAllele(oldAllele);
                if(oldAlleleIndex == -1)
                    throw "missing old allele in likelihood collection";
                if(oldToNewAlleleIndexMap->operator[](oldAlleleIndex) != -1)
                    throw "collision: two new alleles make reference to the same old allele";
                oldToNewAlleleIndexMap->operator[](oldAlleleIndex) = newIndex;
            }
         }
         return oldToNewAlleleIndexMap;
     }

     // pay attention! The data type of parameters is different from the methods above
     shared_ptr<vector<int>> oldToNewAlleleIndexMap(const shared_ptr<std::map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Allele>>>>>& newToOldAlleleMap, int oldAlleleCount,  vector<shared_ptr<Allele>>& newAlleles)
     {
         for(auto& h : newAlleles)
             assert(h);

         auto oldToNewAlleleIndexMap = make_shared<vector<int>>(oldAlleleCount, -1); // -1 indicate that there is no new allele that make reference to that old one.
         for(int newIndex = 0; newIndex < newAlleles.size(); newIndex++)
         {
             auto newAllele = newAlleles[newIndex];
             for(shared_ptr<Allele> oldAllele : *newToOldAlleleMap->at(newAllele))
             {
                 int oldAlleleIndex = indexOfAllele(oldAllele);
                 if(oldAlleleIndex == -1)
                     throw "missing old allele in likelihood collection";
                 if(oldToNewAlleleIndexMap->operator[](oldAlleleIndex) != -1)
                     throw "collision: two new alleles make reference to the same old allele";
                 oldToNewAlleleIndexMap->operator[](oldAlleleIndex) = newIndex;
             }
         }
         return oldToNewAlleleIndexMap;
     }

     /**
   * Returns the index of an allele within the likelihood collection.
   *
   * @param allele the query allele.
   *
   * @throws IllegalArgumentException if {@code allele} is {@code null}.
   *
   * @return -1 if the allele is not included, 0 or greater otherwise.
   */
    int indexOfAllele(shared_ptr<Allele>& allele)
    {
        for(int i=0; i<alleles->size(); i++)
        {
            if(allele.get() == alleles->operator[](i).get())
                return i;
        }
        return -1;
    }

    int indexOfAllele(shared_ptr<Allele>&& allele)
    {
        for(int i=0; i<alleles->size(); i++)
        {
            if(allele.get() == alleles->operator[](i).get())
                return i;
        }
        return -1;
    }

    shared_ptr<vector<vector<int>>> overlappingEvidenceIndicesBySampleIndex(shared_ptr<SimpleInterval>& overlap)
    {
        if(overlap == nullptr)
            return nullptr;

        int sampleCount = samples.size();
        auto result = make_shared<vector<vector<int>>>(sampleCount, vector<int>());


        for(int s=0; s < sampleCount; s++)
        {
            auto& buffer = result->operator[](s);
            buffer.clear();
            vector<shared_ptr<E>>& sampleEvidence = evidenceBySampleIndex->operator[](s);
            int sampleEvidenceCount = sampleEvidence.size();
            buffer.reserve(sampleEvidenceCount);
            for (int r = 0; r < sampleEvidenceCount; r++) {
                if (sampleEvidence[r]->overlaps(overlap)) {
                    buffer.template emplace_back(r);
                }
            }
        }
        return result;
    }

    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    shared_ptr<vector<vector<vector<double>>>> marginalLikelihoods(int oldAlleleCount, int newAlleleCount, vector<int>& oldToNewAlleleIndexMap, const shared_ptr<vector<vector<int>>>& evidenceToKeep)
    {
        int sampleCount = samples.size();
        shared_ptr<vector<vector<vector<double>>>> result = make_shared<vector<vector<vector<double>>>>(sampleCount, vector<vector<double>>());

        for (int s = 0; s < sampleCount; s++) {
            int sampleEvidenceCount = evidenceBySampleIndex->operator[](s).size();
            auto& oldSampleValues = valuesBySampleIndex->operator[](s);
            vector<int>* sampleEvidenceToKeep = evidenceToKeep == nullptr || evidenceToKeep->operator[](s).size() == sampleEvidenceCount ? nullptr : &(evidenceToKeep->operator[](s));
            int newSampleEvidenceCount = sampleEvidenceToKeep == nullptr ? sampleEvidenceCount : sampleEvidenceToKeep->size();
            result->operator[](s) = vector<vector<double>>(newAlleleCount, vector<double>(newSampleEvidenceCount, -std::numeric_limits<double>::infinity()));
            auto& newSampleValues = result->operator[](s);

            // For each old allele and unit of evidence we update the new table keeping the maximum likelihood.
            for (int r = 0; r < newSampleEvidenceCount; r++) {
                for (int a = 0; a < oldAlleleCount; a++) {
                    int oldEvidenceIndex = newSampleEvidenceCount == sampleEvidenceCount ? r : (*sampleEvidenceToKeep)[r];
                    int newAlleleIndex = oldToNewAlleleIndexMap[a];
                    if (newAlleleIndex == -1) {
                        continue;
                    }
                    double likelihood = oldSampleValues[a][oldEvidenceIndex];
                    if (likelihood > newSampleValues[newAlleleIndex][r]) {
                        newSampleValues[newAlleleIndex][r] = likelihood;
                    }
                }
            }
        }
        return result;
    }

    /**
   * Returns the total count of evidence in the evidence-likelihood collection.
   */
   int evidenceCount(){
       int sum = 0;
       for(auto& evidences : *evidenceBySampleIndex)
       {
           sum += evidences.size();
       }
       return sum;
   }
};

template <typename E, typename A>
class SampleMatrix{
private:
    int sampleIndex;
    AlleleLikelihoods<E, A> * likelihood;

public:
    SampleMatrix()
    {
        sampleIndex = 0;
        likelihood = nullptr;
    }

    SampleMatrix(int sampleIndex, AlleleLikelihoods<E, A> * likelihoods)
    {
        this->sampleIndex = sampleIndex;
        this->likelihood = likelihoods;
    }

    void set(int alleleIndex, int evidenceIndex, double value){
        assert(likelihood != nullptr);
        auto& valuesBySampleIndex = *(likelihood->valuesBySampleIndex);
        valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex] = value;
    }

    double get(int alleleIndex, int evidenceIndex){
        assert(likelihood != nullptr);
        auto& valuesBySampleIndex = *(likelihood->valuesBySampleIndex);
        return valuesBySampleIndex[sampleIndex][alleleIndex][evidenceIndex];
    }

    vector<shared_ptr<E>> & evidence(){
        return likelihood->sampleEvidence(sampleIndex);
    }

    vector<shared_ptr<A>>& alleles(){
        return likelihood->getAlleles();
    }

    // different from alleles(), this method returns a shared_ptr
    shared_ptr<vector<shared_ptr<A>>> getAlleles()
    {
        return likelihood->getAlleleList();
    }

    shared_ptr<A> getAllele(int alleleIndex)
    {
        return likelihood->getAllele(alleleIndex);
    }

    int indexOfAllele(shared_ptr<A> allele)
    {
        return likelihood->indexOfAllele(allele);
    }

    AlleleLikelihoods<E, A>* getLikelihoods(){
        assert(likelihood != nullptr);
        return likelihood;
    }

    int numberOfAlleles(){
        return alleles().size();
    }

    int evidenceCount(){
        return likelihood->evidenceBySampleIndex->operator[](sampleIndex).size();
    }

    const vector<vector<double>>& getValuseBySampleIndex()
    {
        return likelihood->valuesBySampleIndex->operator[](sampleIndex);
    }
};

template <typename E, typename A>
class BestAllele{
public:
    /**
     * Null if there is no possible match (no allele?).
     */
    shared_ptr<A> allele;

    /**
     * The containing sample.
     */
    string sample;

    /**
     * The query evidence.
     */
    shared_ptr<E> evidence;

    /**
     * If allele != null, the indicates the likelihood of the evidence.
     */
    double likelihood;

    AlleleLikelihoods<E, A> * alleleLikelihoods;

    /**
     * Confidence that the evidence actually was generated under that likelihood.
     * This is equal to the difference between this and the second best allele match.
     */
    double confidence;

    BestAllele(AlleleLikelihoods<E, A> * likelihoods, int sampleIndex, int evidenceIndex, int bestAlleleIndex,
               double likelihood, double secondBestLikelihood){
        allele = bestAlleleIndex == -1 ? nullptr : likelihoods->alleles->operator[](bestAlleleIndex);
        this->likelihood = likelihood;
        sample = likelihoods->samples[sampleIndex];
        evidence = (*likelihoods->evidenceBySampleIndex)[sampleIndex][evidenceIndex];
        confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
        alleleLikelihoods = likelihoods;
    }

    bool isInformative() {
        return confidence > AlleleLikelihoods<E, A>::LOG_10_INFORMATIVE_THRESHOLD;
    }
};

#endif //MUTECT2CPP_MASTER_ALLELELIKELIHOODS_H
