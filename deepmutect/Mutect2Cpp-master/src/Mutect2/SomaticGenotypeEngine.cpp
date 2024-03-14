//
// Created by lhh on 5/20/22.
//

#include <haplotypecaller/AssemblyBasedCallerUtils.h>
#include <utils/variant/GATKVariantContextUtils.h>
#include "SomaticGenotypeEngine.h"
#include "NaturalLogUtils.h"
#include "SomaticLikelihoodsEngine.h"
#include "variantcontext/builder/GenotypeBuilder.h"
#include "variantcontext/VCFConstants.h"

SomaticGenotypeEngine::SomaticGenotypeEngine(M2ArgumentCollection &MTAC, const string &normalSample,
                                             VariantAnnotatorEngine &annotationEngine, ReferenceCache *cache) : MTAC(
		MTAC), normalSample(normalSample), annotationEngine(annotationEngine), hasNormal(!normalSample.empty()),
                                                                                                                cache(cache) {

}

CalledHaplotypes SomaticGenotypeEngine::callMutations(AlleleLikelihoods<SAMRecord, Haplotype> *logReadLikelihoods,
                                                      AssemblyResultSet &assemblyResultSet,
                                                      ReferenceContext &referenceContext,
                                                      SimpleInterval &activeRegionWindow, SAMFileHeader *header) {

	auto haplotypes = logReadLikelihoods->getAlleles();

	vector<int> startPosKeySet;
	for (int startPosKey: EventMap::buildEventMapsForHaplotypes(haplotypes,
	                                                            assemblyResultSet.getFullReferenceWithPadding(),
	                                                            assemblyResultSet.getFullReferenceWithPaddingLength(),
	                                                            assemblyResultSet.getPaddedReferenceLoc(),
	                                                            MTAC.maxMnpDistance)) {
		if (activeRegionWindow.getStart() <= startPosKey && startPosKey <= activeRegionWindow.getEnd())
			startPosKeySet.push_back(startPosKey);
	}

	auto calledHaplotypes = make_shared<phmap::flat_hash_set<shared_ptr<Haplotype>>>();
	vector<shared_ptr<VariantContext>> returnCalls;

	if (MTAC.likelihoodArgs.phredScaledGlobalReadMismappingRate > 0) {
		logReadLikelihoods->normalizeLikelihoods(
				NaturalLogUtils::qualToLogErrorProb((uint8_t) MTAC.likelihoodArgs.phredScaledGlobalReadMismappingRate));
	}

	AlleleLikelihoods<Fragment, Haplotype> *logFragmentLikelihoods = logReadLikelihoods->groupEvidence(
			&SAMRecord::getName, &Fragment::createAndAvoidFailure);


	for (int loc: startPosKeySet) {
		auto eventsAtThisLoc = AssemblyBasedCallerUtils::getVariantContextsFromActiveHaplotypes(loc, haplotypes, false);
		auto mergedVC = AssemblyBasedCallerUtils::makeMergedVariantContext(eventsAtThisLoc);
		if (mergedVC == nullptr)
			continue;

		auto alleleMapper = AssemblyBasedCallerUtils::createAlleleMapper(mergedVC, loc, haplotypes);
		shared_ptr<AlleleLikelihoods<Fragment, Allele>> logLikelihoods(logFragmentLikelihoods->marginalize(alleleMapper,
		                                                                                                   make_shared<SimpleInterval>(
				                                                                                                   mergedVC->getContig(),
				                                                                                                   mergedVC->getStart(),
				                                                                                                   mergedVC->getEnd())->expandWithinContig(
				                                                                                                   ALLELE_EXTENSION,
				                                                                                                   &header->getSequenceDictionary())));

		vector<SampleMatrix<Fragment, Allele> *> tumorMatrices;
		for (int i = 0; i < logLikelihoods->numberOfSamples(); i++) {
			if (logLikelihoods->getSample(i) != normalSample) {
				tumorMatrices.push_back(logLikelihoods->sampleMatrix(i));
			}
		}
		auto alleleList = tumorMatrices[0];
		SampleMatrix<Fragment, Allele> *logTumorMatrix = combinedLikelihoodMatrix(tumorMatrices, alleleList);
		auto tumorLogOdds = somaticLogOdds(logTumorMatrix);

		vector<SampleMatrix<Fragment, Allele> *> normalMatrices;
		for (int i = 0; i < logLikelihoods->numberOfSamples(); i++) {
			if (logLikelihoods->getSample(i) == normalSample) {
				normalMatrices.push_back(logLikelihoods->sampleMatrix(i));
			}
		}
		SampleMatrix<Fragment, Allele> *logNormalMatrix = combinedLikelihoodMatrix(normalMatrices, alleleList);
		auto normalLogOdds = diploidAltLogOdds(logNormalMatrix);
		auto normalArtifactLogOdds = somaticLogOdds(logNormalMatrix);

		set<shared_ptr<Allele>> forcedAlleles;
		vector<shared_ptr<Allele>> tumorAltAlleles;
		long somaticAltCount = 0;
		for (auto &alternateAllele: mergedVC->getAlternateAlleles()) {
			if (forcedAlleles.find(alternateAllele) != forcedAlleles.end() ||
			    tumorLogOdds->getAlt(alternateAllele) > MTAC.getEmissionLogOdds())
				tumorAltAlleles.emplace_back(alternateAllele);
		}

		for (auto &allele: tumorAltAlleles) {
			if (forcedAlleles.find(allele) != forcedAlleles.end() || !hasNormal || MTAC.genotypeGermlineSites ||
			    normalLogOdds->getAlt(allele) > MathUtils::log10ToLog(MTAC.normalLog10Odds))
				somaticAltCount++;
		}

		// if every alt allele is germline, skip this variant.  However, if some alt alleles are germline and others
		// are not we emit them all so that the filtering engine can see them
		if (somaticAltCount == 0) {
			delete logTumorMatrix->getLikelihoods();
			delete logNormalMatrix->getLikelihoods();
			continue;
		}

		auto allAllelesToEmit = make_shared<vector<shared_ptr<Allele>>>(tumorAltAlleles);
		allAllelesToEmit->emplace_back(mergedVC->getReference());

		vector<shared_ptr<VariantContext>> germlineResourceVariants;
		auto negativeLogPopulationAFAnnotation = getNegativeLogPopulationAFAnnotation(germlineResourceVariants,
		                                                                              tumorAltAlleles,
		                                                                              MTAC.getDefaultAlleleFrequency());
		VariantContextBuilder callVcb(mergedVC);
		callVcb.setAlleles(allAllelesToEmit)->setAttributes(negativeLogPopulationAFAnnotation);

		vector<double> attributeValue;
		for (const auto &a: tumorAltAlleles) {
			attributeValue.push_back(MathUtils::logToLog10(tumorLogOdds->getAlt(a)));
		}
		callVcb.setAttribute(VCFConstants::TUMOR_LOG_10_ODDS_KEY, attributeValue);

		if (hasNormal) {
			callVcb.setAttribute(VCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY,
			                     MathUtils::applyToArrayInPlace(normalArtifactLogOdds->asDoubleArray(tumorAltAlleles),
			                                                    [](double x) { return -MathUtils::logToLog10(x); }));
			callVcb.setAttribute(VCFConstants::NORMAL_LOG_10_ODDS_KEY,
			                     MathUtils::applyToArrayInPlace(normalLogOdds->asDoubleArray(tumorAltAlleles),
			                                                    [](double x) { return MathUtils::logToLog10(x); }));
		}


		addGenotypes(logLikelihoods, allAllelesToEmit, callVcb);
		auto call = callVcb.make();
		auto trimmedCall = GATKVariantContextUtils::trimAlleles(call, true, true);
		auto trimmedAlleles = trimmedCall->getAlleles();
		auto untrimmedAlleles = call->getAlleles();

		auto trimmedToUntrimmedAlleleMap = make_shared<std::map<shared_ptr<Allele>, shared_ptr<vector<shared_ptr<Allele>>>>>();
		for (int i = 0; i < trimmedAlleles.size(); i++) {
			shared_ptr<vector<shared_ptr<Allele>>> untrimmedAllelesList = make_shared<vector<shared_ptr<Allele>>>();
			*untrimmedAllelesList = {untrimmedAlleles[i]};
			trimmedToUntrimmedAlleleMap->insert({trimmedAlleles[i], untrimmedAllelesList});
		}

		auto trimmedLikelihoods = logLikelihoods->marginalize(trimmedToUntrimmedAlleleMap);

		// AlleleLikelihoods for annotation only
		AlleleLikelihoods<SAMRecord, Allele> *logReadAlleleLikelihoods = logReadLikelihoods->marginalize(alleleMapper,
		                                                                                                 make_shared<SimpleInterval>(
				                                                                                                 mergedVC->getContig(),
				                                                                                                 mergedVC->getStart(),
				                                                                                                 mergedVC->getEnd())->expandWithinContig(
				                                                                                                 ALLELE_EXTENSION,
				                                                                                                 &header->getSequenceDictionary()));

		AlleleLikelihoods<SAMRecord, Allele> *trimmedLikelihoodsForAnnotation = logReadAlleleLikelihoods->marginalize(
				trimmedToUntrimmedAlleleMap);

		int len;
		auto refcache = cache->getSubsequenceAt(trimmedCall->getContigInt(), trimmedCall->getStart(),
		                                        trimmedCall->getEnd() + 150, len);
		referenceContext.setCache(refcache, len);
		auto annotatedCall = annotationEngine.annotateContext(trimmedCall, referenceContext,
		                                                      trimmedLikelihoodsForAnnotation);

		for (const auto &allele: call->getAlleles()) {
			for (const auto &h: *alleleMapper->at(allele)) {
				if (h != nullptr)
					calledHaplotypes->insert(h);
			}
		}
		returnCalls.emplace_back(annotatedCall);

		delete logTumorMatrix->getLikelihoods();
		delete logNormalMatrix->getLikelihoods();
		delete trimmedLikelihoods;
		delete logReadAlleleLikelihoods;
		delete trimmedLikelihoodsForAnnotation;
	}

	auto outputCalls = AssemblyBasedCallerUtils::phaseCalls(returnCalls, *calledHaplotypes);
	int eventCount = outputCalls->size();

	auto outputCallsWithEventCountAnnotation = make_shared<std::vector<std::shared_ptr<VariantContext>>>();
	for (auto &vc: *outputCalls) {
		outputCallsWithEventCountAnnotation->emplace_back(
				VariantContextBuilder(vc).setAttribute(VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount)->make());
	}

	delete logFragmentLikelihoods;
	return CalledHaplotypes(outputCallsWithEventCountAnnotation, calledHaplotypes);
}

SampleMatrix<Fragment, Allele> *
SomaticGenotypeEngine::combinedLikelihoodMatrix(vector<SampleMatrix<Fragment, Allele> *> matrices,
                                                SampleMatrix<Fragment, Allele> *alleleList) {
	vector<shared_ptr<Fragment>> reads;
	if (matrices.size() == 1) {
		reads = matrices[0]->evidence();
	} else {
		for (auto matrix: matrices) {
			for (shared_ptr<Fragment> &evidence: matrix->evidence()) {
				reads.emplace_back(evidence);
			}
		}
	}

	vector<string> sample{"COMBINED"};
	map<string, vector<shared_ptr<Fragment>>> evidenceBySample;
	evidenceBySample.insert({"COMBINED", reads});
	auto alleles = alleleList->getAlleles();
	AlleleLikelihoods<Fragment, Allele> *combinedLikelihoods;   // TODO: solve the memory leak problem here
	combinedLikelihoods = new AlleleLikelihoods<Fragment, Allele>(sample, alleles, evidenceBySample);

	int combinedReadIndex = 0;
	auto result = combinedLikelihoods->sampleMatrix(0);
	int alleleCount = result->numberOfAlleles();
	for (auto matrix: matrices) {
		int readCount = matrix->evidenceCount();
		for (int r = 0; r < readCount; r++) {
			for (int a = 0; a < alleleCount; a++) {
				result->set(a, combinedReadIndex, matrix->get(a, r));
			}
			combinedReadIndex++;
		}
	}
	return result;
}

//---this method is slightly different in GATK 4.1.4.1 and GATK 4.2
shared_ptr<PreAlleleCollection<double>>
SomaticGenotypeEngine::somaticLogOdds(SampleMatrix<Fragment, Allele> *logMatrix) {
	int alleleListEnd = logMatrix->alleles().size() - 1;

	auto alleles = logMatrix->alleles();
	for (int i = 0; i < alleleListEnd; i++) {
		if (alleles[i] == Allele::NON_REF_ALLELE)
			throw "<NON_REF> must be last in the allele list.";
	}

	int nonRefIndex = alleles[alleleListEnd] == Allele::NON_REF_ALLELE ? alleleListEnd : -1;
	double logEvidenceWithAllAlleles = logMatrix->evidenceCount() == 0 ? 0.0 : SomaticLikelihoodsEngine::logEvidence(
			logMatrix->getValuseBySampleIndex(), MTAC.minAF, nonRefIndex);

	auto lods = make_shared<PreAlleleCollection<double>>(PreAlleleCollection<double>::Type::ALT_ONLY);
	int refIndex = getRefIndex(logMatrix);

	int numOfAlleles = logMatrix->numberOfAlleles();
	for (int i = 0; i < numOfAlleles; i++) {
		if (i == refIndex)
			continue;

		shared_ptr<Allele> allele = logMatrix->getAllele(i);
		auto logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix<Fragment, Allele>::excludingAllele(logMatrix,
		                                                                                               allele);
		double logEvidenceWithoutThisAllele =
				logMatrixWithoutThisAllele->evidenceCount() == 0 ? 0.0 : SomaticLikelihoodsEngine::logEvidence(
						*getAsRealMatrix(logMatrixWithoutThisAllele), MTAC.minAF,
						logMatrixWithoutThisAllele->numberOfAlleles() > 1 ? nonRefIndex - 1 : -1);
		lods->setAlt(allele, logEvidenceWithAllAlleles - logEvidenceWithoutThisAllele);
	}
	return lods;
}

int SomaticGenotypeEngine::getRefIndex(SampleMatrix<Fragment, Allele> *logMatrix) {
	int numOfAlleles = logMatrix->numberOfAlleles();
	for (int i = 0; i < numOfAlleles; i++) {
		if (logMatrix->getAllele(i)->getIsReference())
			return i;
	}
	throw std::invalid_argument("No ref allele found in likelihoods");
}

int SomaticGenotypeEngine::getRefIndex(SubsettedLikelihoodMatrix<Fragment, Allele> *logMatrix) {
	int numOfAlleles = logMatrix->numberOfAlleles();
	for (int i = 0; i < numOfAlleles; i++) {
		if (logMatrix->getAllele(i)->getIsReference())
			return i;
	}
	throw std::invalid_argument("No ref allele found in likelihoods");
}

shared_ptr<vector<vector<double>>>
SomaticGenotypeEngine::getAsRealMatrix(shared_ptr<SubsettedLikelihoodMatrix<Fragment, Allele>> matrix) {
	assert(matrix != nullptr);
	return getAsRealMatrix(*matrix);
}

shared_ptr<vector<vector<double>>>
SomaticGenotypeEngine::getAsRealMatrix(SubsettedLikelihoodMatrix<Fragment, Allele> &matrix) {
	int numRow = matrix.numberOfAlleles();
	int numColumn = matrix.evidenceCount();
	auto result = make_shared<vector<vector<double>>>(numRow, vector<double>(numColumn));
	for (int i = 0; i < numRow; i++) {
		for (int j = 0; j < numColumn; j++) {
			(*result)[i][j] = matrix.get(i, j);
		}
	}
	return result;
}

shared_ptr<PreAlleleCollection<double>>
SomaticGenotypeEngine::diploidAltLogOdds(SampleMatrix<Fragment, Allele> *matrix) {
	int refIndex = getRefIndex(matrix);
	int numReads = matrix->evidenceCount();
	double homRefLogLikelihood = 0.0;
	for (int i = 0; i < numReads; i++)
		homRefLogLikelihood += matrix->get(refIndex, i);

	auto result = make_shared<PreAlleleCollection<double>>(PreAlleleCollection<double>::Type::ALT_ONLY);
	// hom ref likelihood for the ref allele, het likelihood for alt alleles
	int numAlleles = matrix->numberOfAlleles();
	for (int i = 0; i < numAlleles; i++) {
		if (i == refIndex)
			continue;

		double hetLogLikelihood = 0.0;
		for (int r = 0; r < numReads; r++) {
			vector<double> temp = {matrix->get(refIndex, r), matrix->get(i, r)};
			hetLogLikelihood += NaturalLogUtils::logSumExp(temp) + NaturalLogUtils::LOG_ONE_HALF;
		}
		result->setAlt(matrix->getAllele(i), homRefLogLikelihood - hetLogLikelihood);
	}
	return result;
}

shared_ptr<map<string, vector<double>>>
SomaticGenotypeEngine::getNegativeLogPopulationAFAnnotation(
		vector<shared_ptr<VariantContext>> &germlineResourceVariants,
		vector<shared_ptr<Allele>> &altAlleles,
		double afOfAllelesNotInGermlineResource) {
	shared_ptr<VariantContext> germlineVC = germlineResourceVariants.empty() ? nullptr : germlineResourceVariants[0];
	vector<double> populationAlleleFrequencies = getGermlineAltAlleleFrequencies(altAlleles, germlineVC,
	                                                                             afOfAllelesNotInGermlineResource);
	MathUtils::applyToArrayInPlace(populationAlleleFrequencies, [](double x) {
		return -log10(x);
	});

	auto result = make_shared<map<string, vector<double>>>();
	result->insert({VCFConstants::POPULATION_AF_KEY, populationAlleleFrequencies});
	return result;
}

vector<double> SomaticGenotypeEngine::getGermlineAltAlleleFrequencies(vector<shared_ptr<Allele>> &altAlleles,
                                                                      const shared_ptr<VariantContext> &germlineVC,
                                                                      double afOfAllelesNotInGermlineResource) {
	if (germlineVC != nullptr) {

	}
	return vector<double>(altAlleles.size(), afOfAllelesNotInGermlineResource);
}

void SomaticGenotypeEngine::addGenotypes(const shared_ptr<AlleleLikelihoods<Fragment, Allele>> &logLikelihoods,
                                         const shared_ptr<vector<shared_ptr<Allele>>> &allelesToEmit,
                                         VariantContextBuilder &callVcb) {
	int numberOfSamples = logLikelihoods->numberOfSamples();
	std::vector<std::shared_ptr<Genotype>> genotypes;
	for (int n = 0; n < numberOfSamples; n++) {
		string sample = logLikelihoods->getSample(n);
		SubsettedLikelihoodMatrix logMatrix(logLikelihoods->sampleMatrix(n), allelesToEmit);

		auto alleleCounts = getEffectiveCounts(logMatrix);
		vector<double> flatPriorPseudocounts(logMatrix.numberOfAlleles(), 1.0);
		shared_ptr<vector<double>> alleleFractionsPosterior =
				logMatrix.evidenceCount() == 0 ? make_shared<vector<double>>(flatPriorPseudocounts)
				                               : SomaticLikelihoodsEngine::alleleFractionsPosterior(
						*getAsRealMatrix(logMatrix), flatPriorPseudocounts);
		shared_ptr<vector<double>> tumorAlleleFractionsMean = MathUtils::normalizeSumToOne(alleleFractionsPosterior);

		shared_ptr<Allele> ref = logMatrix.getAllele(getRefIndex(&logMatrix));
		vector<shared_ptr<Allele>> alleles;
		if (sample == normalSample)
			alleles = {ref, make_shared<Allele>(*ref)};
		else
			alleles = *logMatrix.getAlleles();

		// construct AD array for genotypes
		std::vector<int> AD(alleleCounts->size());
		for (int i = 0; i < alleleCounts->size(); i++) {
			AD[i] = (int) round(alleleCounts->operator[](i));
		}

		tumorAlleleFractionsMean->erase((tumorAlleleFractionsMean->end() - 1));
		genotypes.emplace_back(GenotypeBuilder(sample, alleles).setAD(AD).attribute(VCFConstants::ALLELE_FRACTION_KEY,
		                                                                            *tumorAlleleFractionsMean).make());   // TODO: validate it
	}
	callVcb.setGenotypes(genotypes);
}

shared_ptr<vector<double>>
SomaticGenotypeEngine::getEffectiveCounts(SubsettedLikelihoodMatrix<Fragment, Allele> &logLikelihoodMatrix) {
	if (logLikelihoodMatrix.evidenceCount() == 0)
		return make_shared<vector<double>>(logLikelihoodMatrix.numberOfAlleles(), 0.0);

	auto logLikelihoods = getAsRealMatrix(logLikelihoodMatrix);
	return MathUtils::sumArrayFunction(0, logLikelihoods->operator[](0).size(), [&logLikelihoods](int columnIndex) {
		return NaturalLogUtils::normalizeFromLogToLinearSpace(
				SomaticLikelihoodsEngine::getColumnOfLogLikelihood(*logLikelihoods, columnIndex));
	});
}