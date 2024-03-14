//
// Created by 梦想家xixi on 2022/6/6.
//

#ifndef MUTECT2CPP_MASTER_MODEL_H
#define MUTECT2CPP_MASTER_MODEL_H

#include <vector>
#include "samtools/SAMRecord.h"
#include "torch/script.h"
#include "ATen/Parallel.h"
#include "EventMap.h"
#include "AssemblyRegion.h"
#include "variantcontext/VariantContext.h"


class model {
public:
	bool modelRefer(const std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> &reads,
	                std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &allVariantsWithinExtendedRegion,
	                const std::shared_ptr<AssemblyRegion> &regionForGenotyping, ReferenceCache *cache, std::vector<std::string> samplesList, std::string normalSample);

	void Initial(const std::string &modelPath);

	bool isInitialized() const;

private:
	torch::jit::script::Module n_model;

	bool initialized = false;

	static std::vector<std::shared_ptr<SAMRecord>>
	readTrim(std::vector<std::shared_ptr<SAMRecord>> &reads, int start, int end);

	static bool isOverlap(int start, int end, const std::shared_ptr<VariantContext> &vc);

	static std::vector<std::vector<std::vector<int>>>
	generateData(const std::vector<std::shared_ptr<SAMRecord>> &trimCaseReads,
	             const std::vector<std::shared_ptr<SAMRecord>> &trimNormalReads,
	             const std::vector<std::shared_ptr<SAMRecord>> &allReads, int referenceStart,
	             const uint8_t *referenceBases, int referenceBasesLength,
	             int vcStart, int vcEnd, const std::shared_ptr<VariantContext> &vc);

	bool classify(float (*inputs)[6][31]);
};


#endif //MUTECT2CPP_MASTER_MODEL_H
