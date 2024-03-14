//
// Created by cluster on 22-11-5.
//

#include "TumorEvidenceFilter.h"

double TumorEvidenceFilter::calculateErrorProbability(const std::shared_ptr<VariantContext> &vc,
                                                      Mutect2FilteringEngine *filteringEngine,
                                                      std::shared_ptr<ReferenceContext>) {
    std::vector<double> tumorLods = Mutect2FilteringEngine::getTumorLogOdds(vc);
    std::vector<int> ADs = filteringEngine->sumADsOverSamples(vc, true, false);
    int maxIndex = 0;
    int tmp = tumorLods[0];
    for(int i = 0; i < tumorLods.size(); i++) {
        if(tmp < tumorLods[i]) {
            tmp = tumorLods[i];
            maxIndex = i;
        }
    }
    int altCount = ADs[maxIndex+1];
    int totalCount = 0;
    for(int i = 0; i < ADs.size(); i++) {
        totalCount += ADs[i];
    }
    return filteringEngine->getSomaticClusteringModel()->probabilityOfSequencingError(Datum(tumorLods[maxIndex],0, 0, altCount, totalCount, SomaticClusteringModel::indelLength(vc, maxIndex)));
}

std::vector<std::string> TumorEvidenceFilter::requiredAnnotations() {
    return {"TLOD"};
}

ErrorType TumorEvidenceFilter::errorType() {
    return SEQUENCING;
}

std::string TumorEvidenceFilter::filterName() {
    return VCFConstants::TUMOR_EVIDENCE_FILTER_NAME;
}

int TumorEvidenceFilter::filterIndex() {
    return 9;
}
