//
// Created by cluster on 22-11-5.
//

#include "Mutect2FilteringEngine.h"
#include "NaturalLogUtils.h"
#include "MathUtils.h"

#include <utility>
#include "TumorEvidenceFilter.h"
#include "StrandArtifactFilter.h"
#include "FilteredHaplotypeFilter.h"
#include "BaseQualityFilter.h"
#include "MappingQualityFilter.h"
#include "DuplicatedAltReadFilter.h"
#include "PanelOfNormalsFilter.h"
#include "NormalArtifactFilter.h"
#include "NRatioFilter.h"
#include "ReadPositionFilter.h"
#include "MinAlleleFractionFilter.h"
#include "ClusteredEventsFilter.h"
#include "GermlineFilter.h"
#include "MultiallelicFilter.h"
#include "FragmentLengthFilter.h"
#include "PolymeraseSlippageFilter.h"


double Mutect2FilteringEngine::roundFinitePrecisionErrors(double probability) {
    return std::max(std::min(probability, 1.0), 0.0);
}

std::vector<double> Mutect2FilteringEngine::getTumorLogOdds(const std::shared_ptr<VariantContext> &vc) {
    if(vc->hasAttribute(VCFConstants::TUMOR_LOG_10_ODDS_KEY)) {
        auto atrs = vc->getAttributes();
        AttributeValue atr = atrs.at(VCFConstants::TUMOR_LOG_10_ODDS_KEY);
        std::vector<double> result = atr.getAttributeAsDoubleVector();
        for(int i = 0; i < result.size(); i++) {
            result[i] = MathUtils::log10ToLog(result[i]);
        }
        return result;
    } else {
        return {};
    }
}

std::vector<int>
Mutect2FilteringEngine::sumADsOverSamples(const std::shared_ptr<VariantContext> &vc, bool includeTumor, bool includeNormal) {
    std::vector<int> ADs = std::vector<int>(vc->getNAlleles(), 0);
    std::shared_ptr<GenoTypesContext> genotype = vc->getGenotypes();
    std::vector<std::shared_ptr<Genotype>> * genotypes = genotype->getGenotypes();
    for(auto &gt : *genotypes) {
        if((includeTumor) && isTumor(gt.get()) || (includeNormal && isNormal(gt.get()))) {
			vector<int> ads = gt->getAD();
            for(int i = 0; i < vc->getNAlleles(); i++) {
                ADs[i] += ads[i];
            }
        }
    }
    return ADs;
}

bool Mutect2FilteringEngine::isNormal(Genotype *genotype) {
    std::string gt_name = genotype->getSampleName();
    if(normalSample.find(gt_name) != std::string::npos) {
        return true;
    }
    return false;
}

bool Mutect2FilteringEngine::isTumor(Genotype *genotype) {
    return !isNormal(genotype);
}

SomaticClusteringModel* Mutect2FilteringEngine::getSomaticClusteringModel() {
    return &somaticClusteringModel;
}

Mutect2FilteringEngine::Mutect2FilteringEngine(M2FiltersArgumentCollection &MTFAC, std::string  normal) : somaticClusteringModel(MTFAC), normalSample(std::move(normal)),
                                                                                                                thresholdCalculator(MTFAC.initialPosteriorThreshold, MTFAC.maxFalsePositiveRate, MTFAC.fScoreBeta){
    filters.emplace_back(new TumorEvidenceFilter());
    filters.emplace_back(new StrandArtifactFilter());
    filters.emplace_back(new FilteredHaplotypeFilter(100));
    filters.emplace_back(new BaseQualityFilter(MTFAC.minMedianBaseQuality));
    filters.emplace_back(new MappingQualityFilter(MTFAC.minMedianMappingQuality, MTFAC.longIndelLength));
    filters.emplace_back(new DuplicatedAltReadFilter(MTFAC.uniqueAltReadCount));
    filters.emplace_back(new PanelOfNormalsFilter());
    filters.emplace_back(new NormalArtifactFilter());
    filters.emplace_back(new NRatioFilter(MTFAC.nRatio));
    filters.emplace_back(new ReadPositionFilter(MTFAC.minMedianReadPosition));
    filters.emplace_back(new MinAlleleFractionFilter(MTFAC.minAf));
    filters.emplace_back(new ClusteredEventsFilter(2));
    filters.emplace_back(new GermlineFilter());
    filters.emplace_back(new MultiallelicFilter(1));
    filters.emplace_back(new FragmentLengthFilter(10000));
    filters.emplace_back(new PolymeraseSlippageFilter(MTFAC.minSlippageLength, MTFAC.slippageRate));
}

std::vector<int>
Mutect2FilteringEngine::sumStrandCountsOverSamples(const std::shared_ptr<VariantContext> &vc, bool includeTumor,
                                                   bool includeNormal) {
    std::vector<int> result = std::vector<int>(4, 0);
    std::shared_ptr<GenoTypesContext> genotype = vc->getGenotypes();
    std::vector<std::shared_ptr<Genotype>> * genotypes = genotype->getGenotypes();
    for(auto &gt : *genotypes) {
        if((includeTumor) && isTumor(gt.get()) || (includeNormal && isNormal(gt.get()))) {
            if(gt->hasExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY)) {
                std::vector<int> tmp = gt->getExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY).getAttributeAsIntVector();
                for(int i = 0; i < 4; i++) {
                    result[i] += tmp[i];
                }
            }
        }
    }
    return result;
}

Mutect2FilteringEngine::~Mutect2FilteringEngine() {
    for(auto filter : filters) {
        delete filter;
    }
}

void Mutect2FilteringEngine::accumulateData(const std::shared_ptr<VariantContext> &vc,
                                            std::shared_ptr<ReferenceContext> referenceContext) {
    bool flag = false;
    for(auto & allele : vc->getAlleles()) {
        if(allele->getIsNonReference() && !allele->getIsNonRefAllele()) {
            flag = true;
        }
    }
    if(!flag) {
        return;
    }
    ErrorProbabilities errorProbabilities = ErrorProbabilities(filters, vc, this, referenceContext);
    for(auto filter : filters) {
        filter->accumulateDataForLearning(vc, errorProbabilities, this);
    }
    std::vector<int> tumorADs = sumADsOverSamples(vc, true, false);
    std::vector<double> tumorLogOdds = Mutect2FilteringEngine::getTumorLogOdds(vc);
    somaticClusteringModel.record(tumorADs, tumorLogOdds, errorProbabilities.getTechnicalArtifactProbability(), errorProbabilities.getNonSomaticProbability(), vc);
    thresholdCalculator.addArtifactProbability(errorProbabilities.getErrorProbability());
}

void Mutect2FilteringEngine::learnParameters() {
    for(auto & filter : filters) {
        filter->learnParametersAndClearAccumulatedData();
    }
    somaticClusteringModel.learnAndClearAccumulatedData();
    thresholdCalculator.relearnThresholdAndClearAcumulatedProbabilities();
}

double Mutect2FilteringEngine::posteriorProbabilityOfNormalArtifact(double negativeLogOddsOfNormalArtifact) {
    return posteriorProbabilityOfError(negativeLogOddsOfNormalArtifact, somaticClusteringModel.getLogPriorOfVariantVersusArtifact());
}

double Mutect2FilteringEngine::posteriorProbabilityOfError(double logOddsOfRealVersusError, double logPriorOfReal) {
    std::vector<double> unweightedPosteriorOfRealAndError = {logOddsOfRealVersusError + logPriorOfReal,
                                                             NaturalLogUtils::log1mexp(logPriorOfReal)};
    auto posteriorOfRealAndError = NaturalLogUtils::normalizeFromLogToLinearSpace(std::make_shared<std::vector<double>>(unweightedPosteriorOfRealAndError));
    return posteriorOfRealAndError->operator[](1);
}

std::vector<double> Mutect2FilteringEngine::weightedAverageOfTumorAFs(const std::shared_ptr<VariantContext> &vc) {
    double totalWeight = 0;
    std::vector<double> AFs = std::vector<double>(vc->getNAlleles() - 1, 0);
    for(auto & k : *vc->getGenotypes()->getGenotypes()) {
        if(isTumor(k.get())) {
            double weight = 0;
            vector<int> ADs = k->getAD();
            for(int i = 0; i < ADs.size(); i++) {
                weight += ADs[i];
            }
            totalWeight += weight;
            std::vector<double> sampleAFs = {0.0};
            if(k->hasExtendedAttribute(VCFConstants::ALLELE_FREQUENCY_KEY)) {
                sampleAFs = k->getExtendedAttribute(VCFConstants::ALLELE_FREQUENCY_KEY).getAttributeAsDoubleVector();
            }
            for(int i = 0; i < sampleAFs.size(); i++) {
                sampleAFs[i] = sampleAFs[i] * weight;
            }
            for(int i = 0; i < AFs.size(); i++) {
                AFs[i] = sampleAFs[i] + AFs[i];
            }
        }
    }
    for(int i = 0; i < AFs.size(); i++) {
        AFs[i] = AFs[i] / totalWeight;
    }
    return AFs;
}

double Mutect2FilteringEngine::getLogSomaticPrior(const std::shared_ptr<VariantContext> & vc, int altIndex) {
    return somaticClusteringModel.getLogPriorOfSomaticVariant(vc, altIndex);
}

bool
Mutect2FilteringEngine::applyFiltersAndAccumulateOutputStats(const std::shared_ptr<VariantContext> &vc, const std::shared_ptr<ReferenceContext>& referenceContext) {
    bool flag = true;
    ErrorProbabilities errorProbabilities = ErrorProbabilities(filters, vc, this, referenceContext);
    for(auto k : errorProbabilities.getProbabilitiesByFilter()) {
        if(k.second > EPSILON && k.second > thresholdCalculator.getThredshold() - EPSILON) {
            flag = false;
            vc->addFilter(k.first->filterIndex());
        }
    }
    return flag;
}

double Mutect2FilteringEngine::posteriorProbabilityOfError(const std::shared_ptr<VariantContext> &vc, double logOddsOfRealVersusError,
                                                           int altIndex) {
    return posteriorProbabilityOfError(logOddsOfRealVersusError, getLogSomaticPrior(vc, altIndex));
}

