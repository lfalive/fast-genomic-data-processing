//
// Created by cluster on 22-11-9.
//

#ifndef MUTECT2CPP_MASTER_UNIVARIATEPOINTVALUEPAIR_H
#define MUTECT2CPP_MASTER_UNIVARIATEPOINTVALUEPAIR_H


class UnivariatePointValuePair {
private:
    double point;
    double value;

public:
    static UnivariatePointValuePair default_pair;
    UnivariatePointValuePair(double point, double value);
    double getPoint() const {
        return point;
    }
    double getValue() const {
        return value;
    }
};


#endif //MUTECT2CPP_MASTER_UNIVARIATEPOINTVALUEPAIR_H
