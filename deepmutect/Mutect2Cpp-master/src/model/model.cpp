//
// Created by 梦想家xixi on 2022/6/6.
//

#include <locale>
#include "model.h"


std::vector<std::vector<std::vector<int>>>
model::generateData(const std::vector<std::shared_ptr<SAMRecord>> &trimCaseReads,
                    const std::vector<std::shared_ptr<SAMRecord>> &trimNormalReads,
                    const std::vector<std::shared_ptr<SAMRecord>> &allReads, int referenceStart,
                    const uint8_t *referenceBases, int referenceBasesLength, int vcStart, int vcEnd,
                    const std::shared_ptr<VariantContext> &vc) {
	int insertion[31]{0};
	if (referenceStart < 0)
		referenceStart = 0;
	for (const std::shared_ptr<SAMRecord> &read: allReads) {
		std::vector<CigarElement> readCigar = read->getCigarElements();
		int readStart = read->getStart();
		for (CigarElement readCigarElement: readCigar) {
			switch (readCigarElement.getOperator()) {
				case M:
				case D:
				case N:
				case EQ:
				case X:
					readStart += readCigarElement.getLength();
					break;
				case I:
					if (readStart >= vcStart && readStart <= vcEnd) {
						int index = readStart - vcStart;
						if (readCigarElement.getLength() > insertion[index])
							insertion[index] = readCigarElement.getLength();
					}
					break;
				default:
					break;
			}
		}
	}

	int padReferenceLength = 31;
	for (int i: insertion)
		padReferenceLength += i;

	std::vector<std::vector<std::vector<int>>> padOutput = std::vector<std::vector<std::vector<int>>>(5,
	                                                                                                  std::vector<std::vector<int>>(
			                                                                                                  6,
			                                                                                                  std::vector<int>(
					                                                                                                  padReferenceLength,
					                                                                                                  0)));
	int referencePoint = vcStart - referenceStart;
	int referenceOutputPoint = 0;
	if (referenceStart + referenceBasesLength <= vcEnd) {
		for (int i: insertion) {
			if (referenceBasesLength <= referencePoint)
				break;
			if (i == 0) {
				switch (referenceBases[referencePoint]) {
					case 65:
						padOutput[0][0][referenceOutputPoint]++;
						break;
					case 67:
						padOutput[0][1][referenceOutputPoint]++;
						break;
					case 71:
						padOutput[0][2][referenceOutputPoint]++;
						break;
					case 84:
						padOutput[0][3][referenceOutputPoint]++;
						break;
					default:
						padOutput[0][4][referenceOutputPoint]++;
						break;
				}
				referencePoint++;
				referenceOutputPoint++;
			} else {
				for (int k = 0; k < i; k++) {
					padOutput[0][5][referenceOutputPoint + k]++;
				}
				referenceOutputPoint += i;

				switch (referenceBases[referencePoint]) {
					case 65:
						padOutput[0][0][referenceOutputPoint]++;
						break;
					case 67:
						padOutput[0][1][referenceOutputPoint]++;
						break;
					case 71:
						padOutput[0][2][referenceOutputPoint]++;
						break;
					case 84:
						padOutput[0][3][referenceOutputPoint]++;
						break;
					default:
						padOutput[0][4][referenceOutputPoint]++;
						break;
				}
				referencePoint++;
				referenceOutputPoint++;
			}
		}
	} else {
		for (int i: insertion) {
			if (i == 0) {
				switch (referenceBases[referencePoint]) {
					case 65:
						padOutput[0][0][referenceOutputPoint]++;
						break;
					case 67:
						padOutput[0][1][referenceOutputPoint]++;
						break;
					case 71:
						padOutput[0][2][referenceOutputPoint]++;
						break;
					case 84:
						padOutput[0][3][referenceOutputPoint]++;
						break;
					default:
						padOutput[0][4][referenceOutputPoint]++;
						break;
				}
				referencePoint++;
				referenceOutputPoint++;
			} else {
				for (int k = 0; k < i; k++) {
					padOutput[0][5][referenceOutputPoint + k]++;
				}
				referenceOutputPoint += i;

				switch (referenceBases[referencePoint]) {
					case 65:
						padOutput[0][0][referenceOutputPoint]++;
						break;
					case 67:
						padOutput[0][1][referenceOutputPoint]++;
						break;
					case 71:
						padOutput[0][2][referenceOutputPoint]++;
						break;
					case 84:
						padOutput[0][3][referenceOutputPoint]++;
						break;
					default:
						padOutput[0][4][referenceOutputPoint]++;
						break;
				}
				referencePoint++;
				referenceOutputPoint++;
			}
		}
	}

	for (const std::shared_ptr<SAMRecord> &read: trimCaseReads) {
		std::shared_ptr<uint8_t[]> readBases = read->getBases();
		int iteratorInsertion[31]{0};
		std::copy(std::begin(insertion), std::end(insertion), std::begin(iteratorInsertion));
		int readStart = read->getStart();
		int readEnd = read->getEnd();
		int readSoftStart = read->getSoftStart();
		int readBasePoint = readStart - readSoftStart;
		int readsOutputPoint = 0;
		if (vcStart < readStart) {
			for (int i = 0; i < readStart - vcStart; i++)
				readsOutputPoint += iteratorInsertion[i] + 1;
		}
		for (CigarElement readCigarElement: read->getCigarElements()) {
			if (readStart > std::min(readEnd, vcEnd))
				break;
			switch (readCigarElement.getOperator()) {
				case D :
					if (readStart + readCigarElement.getLength() > vcStart) {
						for (int i = std::max(vcStart, readStart);
						     i < std::min(readStart + readCigarElement.getLength(), vcEnd); i++) {
							for (int k = 0; k <= iteratorInsertion[i - vcStart]; k++) {
								padOutput[1][5][readsOutputPoint + k]++;
							}
							readsOutputPoint += iteratorInsertion[i - vcStart] + 1;
						}
					}
					readStart += readCigarElement.getLength();
					break;
				case M:
				case N:
				case EQ:
				case X:
					if (readStart + readCigarElement.getLength() > vcStart) {
						if (vcStart > readStart)
							readBasePoint += vcStart - readStart;
						for (int i = std::max(vcStart, readStart);
						     i < std::min(readStart + readCigarElement.getLength(), vcEnd); i++) {
							for (int k = 0; k < iteratorInsertion[i - vcStart]; k++) {
								padOutput[1][5][readsOutputPoint + k]++;
							}
							readsOutputPoint += iteratorInsertion[i - vcStart];

							switch (readBases[readBasePoint]) {
								case 65:
									padOutput[1][0][readsOutputPoint]++;
									break;
								case 67:
									padOutput[1][1][readsOutputPoint]++;
									break;
								case 71:
									padOutput[1][2][readsOutputPoint]++;
									break;
								case 84:
									padOutput[1][3][readsOutputPoint]++;
									break;
								default:
									padOutput[1][4][readsOutputPoint]++;
									break;
							}
							readBasePoint++;
							readsOutputPoint++;
						}
						readStart += readCigarElement.getLength();
					} else {
						readStart += readCigarElement.getLength();
						readBasePoint += readCigarElement.getLength();
					}
					break;
				case I:
					if (readStart > vcStart) {
						if (iteratorInsertion[readStart - vcStart] <= 0)
							throw std::invalid_argument("out of range");
						for (int i = 0; i < readCigarElement.getLength(); i++) {
							switch (readBases[readBasePoint]) {
								case 65:
									padOutput[1][0][readsOutputPoint]++;
									break;
								case 67:
									padOutput[1][1][readsOutputPoint]++;
									break;
								case 71:
									padOutput[1][2][readsOutputPoint]++;
									break;
								case 84:
									padOutput[1][3][readsOutputPoint]++;
									break;
								default:
									padOutput[1][4][readsOutputPoint]++;
									break;
							}
							readBasePoint++;
							readsOutputPoint++;
						}
						for (int i = readCigarElement.getLength(); i < iteratorInsertion[readStart - vcStart]; i++)
							padOutput[1][5][readsOutputPoint++]++;
						iteratorInsertion[readStart - vcStart] = 0;
					} else {
						readBasePoint += readCigarElement.getLength();
					}
					break;
				default:
					break;
			}
		}
	}

	for (const std::shared_ptr<SAMRecord> &read: trimNormalReads) {
		std::shared_ptr<uint8_t[]> readBases = read->getBases();
		int iteratorInsertion[31]{0};
		std::copy(std::begin(insertion), std::end(insertion), std::begin(iteratorInsertion));
		int readStart = read->getStart();
		int readEnd = read->getEnd();
		int readSoftStart = read->getSoftStart();
		int readBasePoint = readStart - readSoftStart;
		int readsOutputPoint = 0;
		if (vcStart < readStart) {
			for (int i = 0; i < readStart - vcStart; i++)
				readsOutputPoint += iteratorInsertion[i] + 1;
		}
		for (CigarElement readCigarElement: read->getCigarElements()) {
			if (readStart > std::min(readEnd, vcEnd))
				break;
			switch (readCigarElement.getOperator()) {
				case D :
					if (readStart + readCigarElement.getLength() > vcStart) {
						for (int i = std::max(vcStart, readStart);
						     i < std::min(readStart + readCigarElement.getLength(), vcEnd); i++) {
							for (int k = 0; k <= iteratorInsertion[i - vcStart]; k++) {
								padOutput[2][5][readsOutputPoint + k]++;
							}
							readsOutputPoint += iteratorInsertion[i - vcStart] + 1;
						}
					}
					readStart += readCigarElement.getLength();
					break;
				case M:
				case N:
				case EQ:
				case X:
					if (readStart + readCigarElement.getLength() > vcStart) {
						if (vcStart > readStart)
							readBasePoint += vcStart - readStart;
						for (int i = std::max(vcStart, readStart);
						     i < std::min(readStart + readCigarElement.getLength(), vcEnd); i++) {
							for (int k = 0; k < iteratorInsertion[i - vcStart]; k++) {
								padOutput[2][5][readsOutputPoint + k]++;
							}
							readsOutputPoint += iteratorInsertion[i - vcStart];

							switch (readBases[readBasePoint]) {
								case 65:
									padOutput[2][0][readsOutputPoint]++;
									break;
								case 67:
									padOutput[2][1][readsOutputPoint]++;
									break;
								case 71:
									padOutput[2][2][readsOutputPoint]++;
									break;
								case 84:
									padOutput[2][3][readsOutputPoint]++;
									break;
								default:
									padOutput[2][4][readsOutputPoint]++;
									break;
							}
							readBasePoint++;
							readsOutputPoint++;
						}
						readStart += readCigarElement.getLength();
					} else {
						readStart += readCigarElement.getLength();
						readBasePoint += readCigarElement.getLength();
					}
					break;
				case I:
					if (readStart > vcStart) {
						if (iteratorInsertion[readStart - vcStart] <= 0)
							throw std::invalid_argument("out of range");
						for (int i = 0; i < readCigarElement.getLength(); i++) {
							switch (readBases[readBasePoint]) {
								case 65:
									padOutput[2][0][readsOutputPoint]++;
									break;
								case 67:
									padOutput[2][1][readsOutputPoint]++;
									break;
								case 71:
									padOutput[2][2][readsOutputPoint]++;
									break;
								case 84:
									padOutput[2][3][readsOutputPoint]++;
									break;
								default:
									padOutput[2][4][readsOutputPoint]++;
									break;
							}
							readBasePoint++;
							readsOutputPoint++;
						}
						for (int i = readCigarElement.getLength(); i < iteratorInsertion[readStart - vcStart]; i++)
							padOutput[2][5][readsOutputPoint++]++;
						iteratorInsertion[readStart - vcStart] = 0;
					} else {
						readBasePoint += readCigarElement.getLength();
					}
					break;
				default:
					break;
			}
		}
	}

	int matricStart = 0;

	for (int i = 0; i < 16; i++) {
		matricStart += insertion[i];
	}

	std::vector<std::vector<std::vector<int>>> result = std::vector<std::vector<std::vector<int>>>(4,
	                                                                                               std::vector<std::vector<int>>(
			                                                                                               6,
			                                                                                               std::vector<int>(
					                                                                                               31,
					                                                                                               0)));
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 6; j++)
			for (int k = 0; k < 31; k++) {
				result[i][j][k] = padOutput[i][j][matricStart + k];
			}

	for (int i = 0; i < 31; i++)
		result[3][0][i] = insertion[i];
	return result;
}

std::vector<std::shared_ptr<SAMRecord>>
model::readTrim(std::vector<std::shared_ptr<SAMRecord>> &reads, int start, int end) {
	std::vector<std::shared_ptr<SAMRecord>> trimReads;
	for (const std::shared_ptr<SAMRecord> &read: reads) {
		if (read->getEnd() < start || read->getStart() > end)
			continue;
		trimReads.emplace_back(read);
	}
	return trimReads;
}

bool model::isOverlap(int start, int end, const std::shared_ptr<VariantContext> &vc) {
	return vc->getStart() >= start && std::max(vc->getEnd(), vc->getStart() + 8) <= end;
}

bool model::modelRefer(const std::shared_ptr<std::map<std::string, std::vector<std::shared_ptr<SAMRecord>>>> &reads,
                       std::set<std::shared_ptr<VariantContext>, VariantContextComparator> &allVariantsWithinExtendedRegion,
                       const std::shared_ptr<AssemblyRegion> &regionForGenotyping, ReferenceCache *cache, std::vector<std::string> samplesList, std::string normalSample) {
	if (allVariantsWithinExtendedRegion.empty()) // no variants,
		return false;

	std::vector <std::shared_ptr<VariantContext>> withinActiveRegion;
	std::shared_ptr <SimpleInterval> originalRegionRange = regionForGenotyping->getSpan();

	std::vector <std::shared_ptr<SAMRecord>> caseReads, normalReads;
	for (const auto &item: samplesList) {
		if (item != normalSample)
			caseReads = reads->at(item);
		else
			normalReads = reads->at(normalSample);
	}
	int referenceStart = regionForGenotyping->getExtendedSpan()->getStart() - 15;
	int ref_len = 0;
	std::shared_ptr < uint8_t[] > referenceBases = regionForGenotyping->getFullReference(cache, 15, ref_len);
	for (const std::shared_ptr <VariantContext> &vc: allVariantsWithinExtendedRegion) {
		std::shared_ptr <SimpleInterval> vcLoc = std::make_shared<SimpleInterval>(vc->getContig(), vc->getStart(),
		                                                                          vc->getEnd());
		if (originalRegionRange->overlaps(vcLoc))
			withinActiveRegion.emplace_back(vc);
	}
	int position = 0;
	for (const std::shared_ptr <VariantContext> &vc: withinActiveRegion) {
		if (position > vc->getEnd())
			continue;
		int vcStart = vc->getStart() - 15;
		int vcEnd = vc->getStart() + 15;
		if (vcStart < 0) {
			vcStart = 0;
			vcEnd = 30;
		}
		std::vector <std::shared_ptr<SAMRecord>> trimCaseReads = readTrim(caseReads, vcStart, vcEnd);
		std::vector <std::shared_ptr<SAMRecord>> trimNormalReads = readTrim(normalReads, vcStart, vcEnd);
		std::vector <std::shared_ptr<SAMRecord>> allReads;
		for (const auto &read: caseReads) {
			allReads.emplace_back(read);
		}
		for (const auto &read: normalReads) {
			allReads.emplace_back(read);
		}
		std::vector < std::vector < std::vector < int>>> result = generateData(trimCaseReads, trimNormalReads, allReads,
		                                                                       referenceStart, referenceBases.get(),
		                                                                       ref_len,
		                                                                       vcStart, vcEnd, vc);
		int count2 = 15;
		int index2 = 15;
		while (count2 < 30) {
			count2 += result[3][0][index2 + 1] + 1;
			index2++;
		}
		if (index2 == 29)
			index2++;
		int regionEnd = vcStart + index2;
		position = regionEnd;
		int count[2][31]{0};
		float inputs[5][6][31]{0};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 31; j++) {
				count[i][j] = result[i + 1][0][j] + result[i + 1][1][j] + result[i + 1][2][j] + result[i + 1][3][j] +
				              result[i + 1][4][j] + result[i + 1][5][j];
			}
		}
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 31; j++) {
				inputs[0][i][j] = static_cast<float>(result[0][i][j]);
			}
		}
		for (int i = 1; i < 3; i++) {
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 31; k++) {
					inputs[i][j][k] = static_cast<float>(result[i][j][k]) /
					                  (static_cast<float>(count[i - 1][k]) + 0.00000000001f);
				}
			}
		}
		for (int i = 3; i < 5; i++) {
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 31; k++) {
					inputs[i][j][k] = static_cast<float>(result[i - 2][j][k]) /
					                  (static_cast<float>(result[1][j][k]) + static_cast<float>(result[2][j][k]
					                                                                            + 0.00000000001f));
				}
			}
		}

		if (classify(inputs))
			return true;
	}
	return false;
}

void model::Initial(const std::string &modelPath) {
	try {
		n_model = torch::jit::load(modelPath);
		at::set_num_threads(1);
		initialized = true;
	}
	catch (const c10::Error &e) {
		std::cerr << "error loading the model\n";
		initialized = false;
	}
}

bool model::classify(float (*fre)[6][31]) {
	try {
		n_model.eval();
		torch::Tensor inputs = torch::from_blob(fre, {1, 30, 31}, torch::kFloat);
		inputs = inputs.transpose(1, 2);
		torch::Tensor out_tensor = n_model.forward({inputs}).toTensor();
		bool output = out_tensor[0][0].item().toFloat() > 0.9999999995;
		//std::cout << output << std::endl;
		return output;
	}
	catch (const c10::Error &e) {
		std::cerr << "error loading the model\n";
		initialized = false;
		return false;
	}
}

bool model::isInitialized() const {
	return initialized;
}
