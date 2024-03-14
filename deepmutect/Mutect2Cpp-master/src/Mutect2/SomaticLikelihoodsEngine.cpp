//
// Created by lhh on 5/28/22.
//

#include <boost/math/special_functions/beta.hpp>
#include "SomaticLikelihoodsEngine.h"
#include "Mutect2Utils.h"
#include "utils/Dirichlet.h"
#include "MathUtils.h"
#include "NaturalLogUtils.h"

double SomaticLikelihoodsEngine::logEvidence(const std::vector<std::vector<double>> &logLikelihoods, double minAF,
                                             int nonRefIndex) {
	int RowDimension = logLikelihoods.size();
	std::vector<double> flatPrior(RowDimension, 1.0);
	return logEvidence(logLikelihoods, flatPrior, minAF, nonRefIndex);
}

double SomaticLikelihoodsEngine::logEvidence(const std::vector<std::vector<double>> &logLikelihoods,
                                             std::vector<double> &priorPseudocounts, double alleleFractionThreshold,
                                             int nonRefIndex) {
	int numberOfAlleles = logLikelihoods.size();
	if (numberOfAlleles != priorPseudocounts.size())
		throw std::invalid_argument("Must have one pseudocount per allele.");

	auto alleleFractionsPosteriorArray = alleleFractionsPosterior(logLikelihoods, priorPseudocounts);
	double priorContribution = logDirichletNormalization(priorPseudocounts);
	double posteriorContribution = -logDirichletNormalization(*alleleFractionsPosteriorArray);
	double posteriorTotal = MathUtils::sum(alleleFractionsPosteriorArray);
	double thresholdedPosteriorContribution = posteriorContribution;
	if (nonRefIndex > 0) {
		thresholdedPosteriorContribution += log(1 - boost::math::ibeta(alleleFractionThreshold,
		                                                               alleleFractionsPosteriorArray->operator[](
				                                                               nonRefIndex), posteriorTotal -
		                                                                                     alleleFractionsPosteriorArray->operator[](
				                                                                                     nonRefIndex)));
	}

	auto logAlleleFractions = Dirichlet(alleleFractionsPosteriorArray).effectiveLogMultinomialWeights();
	double likelihoodsAndEntropyContribution = 0.0;
	int columnNums = logLikelihoods[0].size();
	for (int i = 0; i < columnNums; i++) {
		auto logLikelihoodsForRead = getColumnOfLogLikelihood(logLikelihoods, i);
		auto responsibilities = NaturalLogUtils::posteriors(*logAlleleFractions, *logLikelihoodsForRead);
		double entropyContribution = MathUtils::sum(
				MathUtils::applyToArray(responsibilities, SomaticLikelihoodsEngine::xLogx));
		likelihoodsAndEntropyContribution +=
				likelihoodsContribution(logLikelihoodsForRead, responsibilities) - entropyContribution;
	}

	return priorContribution + thresholdedPosteriorContribution + likelihoodsAndEntropyContribution;
}

std::shared_ptr<std::vector<double>>
SomaticLikelihoodsEngine::alleleFractionsPosterior(const std::vector<std::vector<double>> &logLikelihoods,
                                                   std::vector<double> &priorPseudocounts) {
	int numberOfAlleles = logLikelihoods.size();
	if (numberOfAlleles != priorPseudocounts.size())
		throw std::invalid_argument("Must have one pseudocount per allele.");
	auto dirichletPosterior = make_shared<std::vector<double>>(numberOfAlleles, 1.0);
	bool converged = false;

	while (!converged) {
		// alleleCounts = \sum_r \bar{z}_r, where \bar{z}_r is an a-dimensional vector of the expectation of z_r with respect to q(f)
		auto alleleCounts = getEffectiveCounts(logLikelihoods, dirichletPosterior);

		int size = alleleCounts->size();
		auto newDirichletPosterior = make_shared<std::vector<double>>(size);
		for (int i = 0; i < size; i++) {
			newDirichletPosterior->operator[](i) = alleleCounts->operator[](i) + priorPseudocounts[i];
		}
		converged = MathUtils::distance1(*dirichletPosterior, *newDirichletPosterior) < CONVERGENCE_THRESHOLD;
		dirichletPosterior = newDirichletPosterior;
	}
	return dirichletPosterior;
}

std::shared_ptr<std::vector<double>>
SomaticLikelihoodsEngine::getEffectiveCounts(const std::vector<std::vector<double>> &logLikelihoods,
                                             std::shared_ptr<std::vector<double>> dirichletPrior) {
	auto effectiveLogWeights = Dirichlet(dirichletPrior).effectiveLogMultinomialWeights();
	int min = 0, max = logLikelihoods[0].size();
	auto result = NaturalLogUtils::posteriors(effectiveLogWeights, getColumnOfLogLikelihood(logLikelihoods, min));
	for (int n = min + 1; n < max; n++) {
		auto newValues = NaturalLogUtils::posteriors(effectiveLogWeights, getColumnOfLogLikelihood(logLikelihoods, n));
		if (newValues->size() != result->size())
			throw std::invalid_argument("array function returns different sizes for different inputs!");
		for (int i = 0; i < result->size(); i++) {
			result->operator[](i) += newValues->operator[](i);
		}
	}
	return result;
}

std::shared_ptr<std::vector<double>>
SomaticLikelihoodsEngine::getColumnOfLogLikelihood(const std::vector<std::vector<double>> &logLikelihoods,
                                                   int colomnIndex) {
	auto logLikelihoodsCoulmn = make_shared<std::vector<double>>(logLikelihoods.size());
	int numRow = logLikelihoods.size();
	for (int i = 0; i < numRow; i++) {
		logLikelihoodsCoulmn->operator[](i) = logLikelihoods[i][colomnIndex];
	}
	return logLikelihoodsCoulmn;
}

double SomaticLikelihoodsEngine::logDirichletNormalization(std::vector<double> &dirichletParams) {
	double logNumerator = lgamma(MathUtils::sum(dirichletParams));
	double logDenominator = MathUtils::sum(
			MathUtils::applyToArray(dirichletParams, [](double a) { return lgamma(a); }));
	return logNumerator - logDenominator;
}

double SomaticLikelihoodsEngine::xLogx(double x) {
	return x < 1e-8 ? 0.0 : x * log(x);
}

double SomaticLikelihoodsEngine::likelihoodsContribution(const std::shared_ptr<std::vector<double>>& logLikelihoodsForRead,
                                                         const std::shared_ptr<std::vector<double>>& responsibilities) {
	// this is a safe version of MathUtils.sum(MathArrays.ebeMultiply(logLikelihoodsForRead, responsibilities))
	// in case the likelihood is zero, and the log likelihood in -Infinity, we have the responsibility is zero and the
	// contribution x * log(y), where x and y go to zero at the same rate (ie within a constant factor of each other
	// since the responsibility is related to the likelihood via the prior), is undefined but should be treated as zero.
	double result = 0.0;
	for (int n = 0; n < logLikelihoodsForRead->size(); n++) {
		result += (responsibilities->operator[](n) < NEGLIGIBLE_RESPONSIBILITY ? 0 :
		           logLikelihoodsForRead->operator[](n) * responsibilities->operator[](n));
	}
	return result;
}