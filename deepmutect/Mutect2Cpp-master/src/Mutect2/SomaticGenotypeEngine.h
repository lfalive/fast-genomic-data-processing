//
// Created by lhh on 5/20/22.
//

#ifndef MUTECT2CPP_MASTER_SOMATICGENOTYPEENGINE_H
#define MUTECT2CPP_MASTER_SOMATICGENOTYPEENGINE_H

#include "M2ArgumentCollection.h"
#include "VariantAnnotatorEngine.h"
#include "haplotypecaller/CalledHaplotypes.h"
#include "engine/ReferenceContext.h"
#include "PreAlleleCollection.h"
#include "SubsettedLikelihoodMatrix.h"
#include "variantcontext/builder/VariantContextBuilder.h"

class SomaticGenotypeEngine {
private:
    M2ArgumentCollection& MTAC;
    string normalSample;
    bool hasNormal;
    ReferenceCache * cache;

    static shared_ptr<map<string, vector<double>>> getNegativeLogPopulationAFAnnotation(vector<shared_ptr<VariantContext>>& germlineResourceVariants, vector<shared_ptr<Allele>>& altAlleles, double afOfAllelesNotInGermlineResource);

    static vector<double> getGermlineAltAlleleFrequencies(vector<shared_ptr<Allele>>& altAlleles, const shared_ptr<VariantContext>& germlineVC,  double afOfAllelesNotInGermlineResource);

    static shared_ptr<vector<double>> getEffectiveCounts(SubsettedLikelihoodMatrix<Fragment, Allele>& logLikelihoodMatrix);
protected:
    VariantAnnotatorEngine& annotationEngine;

public:
    const static int ALLELE_EXTENSION = 2;

    SomaticGenotypeEngine(M2ArgumentCollection& MTAC, const string& normalSample, VariantAnnotatorEngine& annotationEngine, ReferenceCache * cache);

    /**
    * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
    * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
    *
    * The list of samples we're working with is obtained from the readLikelihoods
    * @param logReadLikelihoods                       Map from reads->(haplotypes,likelihoods)
    * @param activeRegionWindow                     Active window
    * @param withBamOut                            whether to annotate reads in readLikelihoods for future writing to bamout
    * @param emitRefConf                           generate reference confidence (GVCF) data?
    * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
    */
    CalledHaplotypes callMutations(AlleleLikelihoods<SAMRecord, Haplotype>* logReadLikelihoods, AssemblyResultSet& assemblyResultSet, ReferenceContext& referenceContext, SimpleInterval& activeRegionWindow, SAMFileHeader * header);

    static SampleMatrix<Fragment, Allele>* combinedLikelihoodMatrix(vector<SampleMatrix<Fragment, Allele>*> matrices, SampleMatrix<Fragment, Allele>* alleleList);

    // compute the likelihoods that each allele is contained at some allele fraction in the sample
    shared_ptr<PreAlleleCollection<double>> somaticLogOdds(SampleMatrix<Fragment, Allele>* logMatrix);

    static int getRefIndex(SampleMatrix<Fragment, Allele>* logMatrix);

    static int getRefIndex(SubsettedLikelihoodMatrix<Fragment, Allele>* logMatrix);

    static shared_ptr<vector<vector<double>>> getAsRealMatrix(shared_ptr<SubsettedLikelihoodMatrix<Fragment, Allele>> matrix);

    static shared_ptr<vector<vector<double>>> getAsRealMatrix(SubsettedLikelihoodMatrix<Fragment, Allele>& matrix);

    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * @param matrix a matrix of log likelihoods
     */
    shared_ptr<PreAlleleCollection<double>> diploidAltLogOdds(SampleMatrix<Fragment, Allele>* matrix);

    void addGenotypes(const shared_ptr<AlleleLikelihoods<Fragment, Allele>>& logLikelihoods, const shared_ptr<vector<shared_ptr<Allele>>>& allelesToEmit, VariantContextBuilder& callVcb);

    void setReferenceCache(ReferenceCache * cache) {
        this->cache = cache;
    }
};


#endif //MUTECT2CPP_MASTER_SOMATICGENOTYPEENGINE_H
