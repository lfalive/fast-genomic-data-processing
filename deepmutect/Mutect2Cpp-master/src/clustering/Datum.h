//
// Created by cluster on 22-11-7.
//

#ifndef MUTECT2CPP_MASTER_DATUM_H
#define MUTECT2CPP_MASTER_DATUM_H


class Datum {
private:
    double tumorLogOdds;
    double artifactProb;
    double nonSequencingErrorProb;
    int altCount;
    int totalCount;
    int indelLength;

public:
    Datum(double tumorLogOdds, double artifactProb, double nonSomaticProb, int altCount, int totalCount, int indelLength);

    double getTumorLogOdds() const { return tumorLogOdds; }

    double getArtifactProb() const { return artifactProb; }

    double getNonSequencingErrorProb() const { return nonSequencingErrorProb; }

    int getAltCount() const { return altCount; }

    int getTotalCount() const { return totalCount; }

    int getIndelLength() const { return indelLength; }
};


#endif //MUTECT2CPP_MASTER_DATUM_H
