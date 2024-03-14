//
// Created by cluster on 22-11-14.
//

#include "ClusteredEventsFilter.h"

ClusteredEventsFilter::ClusteredEventsFilter(int maxEventsInRegion) : maxEventsInRegion(maxEventsInRegion){

}

ErrorType ClusteredEventsFilter::errorType() {
    return ARTIFACT;
}

bool
ClusteredEventsFilter::isArtifact(const std::shared_ptr<VariantContext> &vc, Mutect2FilteringEngine *filteringEngine) {
    int eventCount = vc->getAttributeAsInt(VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY, -1);
    return eventCount > maxEventsInRegion;
}

std::string ClusteredEventsFilter::filterName() {
    return VCFConstants::CLUSTERED_EVENTS_FILTER_NAME;
}

std::vector<std::string> ClusteredEventsFilter::requiredAnnotations() {
    return {VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY};
}

int ClusteredEventsFilter::filterIndex() {
    return 16;
}
