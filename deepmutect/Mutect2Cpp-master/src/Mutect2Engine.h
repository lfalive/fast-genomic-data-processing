//
// The core algorithm of Mutect2 is here
// Created by lhh on 10/19/21.
//

#ifndef MUTECT2CPP_MASTER_MUTECT2ENGINE_H
#define MUTECT2CPP_MASTER_MUTECT2ENGINE_H

#include <vector>
#include "ActivityProfileState.h"
#include "M2ArgumentCollection.h"
#include "ReferenceCache.h"
#include "engine/ReferenceContext.h"
#include "utils/PeUtils.h"
#include "samtools/SAMFileHeader.h"
#include "engine/AlignmentContext.h"
#include "utils/ReadPileup.h"
#include "ReadCache.h"
#include "model/model.h"
#include "variantcontext/VariantContext.h"
#include "ReadThreadingAssembler.h"
#include "haplotypecaller/MutectReadThreadingAssemblerArgumentCollection.h"
#include "haplotypecaller/AssemblyRegionTrimmer.h"
#include "haplotypecaller/PairHMMLikelihoodCalculationEngine.h"
#include "SmithWatermanAligner.h"
#include "SomaticGenotypeEngine.h"


class Mutect2Engine {
private:
	int minCallableDepth;
	std::vector<string> samplesList;
	std::string &normalSample;
	SAMFileHeader *header;
	M2ArgumentCollection &MTAC;
	ReadThreadingAssembler assemblyEngine;
	MutectReadThreadingAssemblerArgumentCollection assemblerArgs;
	PairHMMLikelihoodCalculationEngine *likelihoodCalculationEngine;
	AssemblyRegionTrimmer trimmer;
	SmithWatermanAligner *aligner;
	SomaticGenotypeEngine genotypingEngine;
	model mymodel;

	static std::shared_ptr<std::vector<char>> altQuals(ReadPileup &pileup, char refBase, int pcrErrorQual);

	static int getCurrentOrFollowingIndelLength(PeUtils &pe);

	static char indelQual(int indelLength);

	static bool isNextToUsefulSoftClip(PeUtils &pe);

	/**
	 * this implement the isActive() algorithm described in docs/mutect/mutect.pdf
	 * the multiplicative factor is for the special case where we pass a singleton list
	 * of alt quals and want to duplicate that alt qual over multiple reads
	 */
	static double logLikelihoodRatio(int refCount, const std::shared_ptr<std::vector<char>> &altQuals);

	static double logLikelihoodRatio(int nRef, const std::shared_ptr<std::vector<char>> &altQuals, int repeatFactor);

	bool hasNormal();

	static void removeUnmarkedDuplicates(const std::shared_ptr<AssemblyRegion> &assemblyRegion);

	static void removeReadStubs(const std::shared_ptr<AssemblyRegion> &assemblyRegion);

public:
	const static int READ_QUALITY_FILTER_THRESHOLD = 20;
	const static int MIN_READ_LENGTH = 30;
	const static int MINIMUM_BASE_QUALITY = 6;
	const static int HUGE_FRAGMENT_LENGTH = 1000000;

	int callableSites;  // in GATK4, this variable is a MutableInt class object
	ReferenceCache *refCache;

	Mutect2Engine(M2ArgumentCollection &MTAC, SAMFileHeader *samFileHeader, const std::string &modelPath,
	              VariantAnnotatorEngine &annotatorEngine);

	~Mutect2Engine();

	std::shared_ptr<ActivityProfileState> isActive(AlignmentContext &context, int refName);

	static void fillNextAssemblyRegionWithReads(const std::shared_ptr<AssemblyRegion> &region, ReadCache &readCache) ;

	std::vector<std::shared_ptr<VariantContext>>
	callRegion(const std::shared_ptr<AssemblyRegion> &originalAssemblyRegion, ReferenceContext &referenceContext);

	// Maybe this variable can be removed in the multi-thread mode
	void setReferenceCache(ReferenceCache *cache);

	static void printVariationContexts(const std::shared_ptr<AssemblyRegion> &region,
	                            const std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &vcs);

	static void printVariationContexts(const std::shared_ptr<AssemblyRegion> &region,
	                            const std::vector<std::shared_ptr<VariantContext>> &vcs);

	static void printVariationContext(const std::shared_ptr<VariantContext> &vc);

	static void
	printReadsMap(const std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> &reads);

protected:
	std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>>
	splitReadsBySample(const std::vector<std::shared_ptr<SAMRecord>> &reads);

};


#endif //MUTECT2CPP_MASTER_MUTECT2ENGINE_H
