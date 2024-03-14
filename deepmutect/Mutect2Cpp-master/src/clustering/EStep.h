//
// Created by cluster on 22-11-8.
//

#ifndef MUTECT2CPP_MASTER_ESTEP_H
#define MUTECT2CPP_MASTER_ESTEP_H


class EStep {
private:
    double forwardArtifactResponsibility;
    double reverseArtifactResponsibility;
    int forwardCount;
    int reverseCount;
    int forwardAltCount;
    int reverseAltCount;

public:
    EStep(double forwardArtifactResponsibility, double reverseArtifactResponsibility, int forwardCount, int reverseCount, int forwardAltCount, int reverseAltCount);
    double getForwardArtifactResponsibility() const {
        return forwardArtifactResponsibility;
    }

    double getReverseArtifactResponsibility() const {
        return reverseArtifactResponsibility;
    }

    double getArtifactProbability() const {
        return getForwardArtifactResponsibility() + getReverseArtifactResponsibility();
    }

    int getForwardCount() const {
        return forwardCount;
    }

    int getReverseCount() const {
        return reverseCount;
    }

    int getForwardAltCount() const {
        return forwardAltCount;
    }

    int getReverseAltCount() const {
        return reverseAltCount;
    }
};


#endif //MUTECT2CPP_MASTER_ESTEP_H
