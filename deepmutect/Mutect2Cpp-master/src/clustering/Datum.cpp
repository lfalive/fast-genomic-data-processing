//
// Created by cluster on 22-11-7.
//

#include "Datum.h"

Datum::Datum(double tumorLogOdds, double artifactProb, double nonSomaticProb, int altCount, int totalCount,
             int indelLength) : tumorLogOdds(tumorLogOdds), artifactProb(artifactProb), altCount(altCount), totalCount(totalCount), indelLength(indelLength){
    nonSequencingErrorProb = 1 - (1 - artifactProb) * (1 - nonSomaticProb);
}
