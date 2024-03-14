/**
 * The implementation of QualQuantizer class
 */

#include <iostream>
#include <math.h>
#include "QualQuantizer.h"
#include "QualityUtils.h"
using namespace std;



QualQuantizer::QualInterval::QualInterval(int qStart, int qEnd, long nObservation, long nErrors, int level, int fixedQual, int minQual)
{
    this->qStart = qStart;
    this->qEnd = qEnd;
    this->nObservations = nObservation;
    this->nErrors = nErrors;
    this->fixedQual = fixedQual;
    this->level = level;
    this->mergeOrder = 0;
    this->subIntervals = nullptr;

    this->minQual = minQual;
}

QualQuantizer::QualInterval::QualInterval(int qStart, int qEnd, long nObservation, long nErrors, int level, int fixedQual, int minQual ,set<QualInterval> * subIntervals)
{
    this->qStart = qStart;
    this->qEnd = qEnd;
    this->nObservations = nObservation;
    this->nErrors = nErrors;
    this->fixedQual = fixedQual;
    this->level = level;
    this->mergeOrder = 0;
    this->subIntervals = subIntervals;

    this->minQual = minQual;
}

QualQuantizer::QualInterval::QualInterval(const QualInterval & interval)
{
    this->qStart = interval.qStart;
    this->qEnd = interval.qEnd;
    this->nObservations = interval.nObservations;
    this->nErrors = interval.nErrors;
    this->fixedQual = interval.fixedQual;
    this->level = interval.level;
    this->mergeOrder = interval.mergeOrder;
    this->subIntervals = interval.subIntervals; //---copy of pointer

    this->minQual = interval.minQual;
    //cout << "copy constructor" << endl;
}

QualQuantizer::QualInterval::~QualInterval()
{
    //cout << "destructor called" << " " << qStart << " - " << qEnd << endl;
    //delete subIntervals;  //---the contents of subIntervals are objects
}

int QualQuantizer::QualInterval::compareTo(QualInterval *qualInterval)
{
    if (this->qStart < qualInterval->qStart)
        return -1;
    else if (this->qStart > qualInterval->qStart)
        return 1;
    else
        return 0;
}

double QualQuantizer::QualInterval::getPenalty()
{
    //return calcPenalty(getErrorRate());
    double Penalty = calcPenalty(getErrorRate());
    return Penalty;
}

double QualQuantizer::QualInterval::calcPenalty(double globalErrorRate) const
{
    if (globalErrorRate == 0.0)
        return 0.0;

    if (this->subIntervals == nullptr || this->subIntervals->empty())
    {
        if (this->qEnd <= this->minQual)
            return 0;
        else
            return abs(log10(getErrorRate()) - log10(globalErrorRate)) * nObservations;
    }else{
        double sum = 0;
        set<QualInterval >::iterator itr;
        for (itr = this->subIntervals->begin(); itr != this->subIntervals->end(); itr++)
        {
            sum += itr->calcPenalty(globalErrorRate);
        }
        return sum;
    }
}

double QualQuantizer::QualInterval::getErrorRate() const
{
    if (fixedQual != -1)
        return QualityUtils::qualToErrorProb((uint8_t)fixedQual);
    else if (nObservations == 0)
        return 0.0;
    else
        return (nErrors + 1) / (1.0 * (nObservations+1));
}

char QualQuantizer::QualInterval::getQual() const
{
    if (fixedQual == -1)
        return QualityUtils::errorProbToQual(getErrorRate());
    else
        return (char)fixedQual;
}

void QualQuantizer::QualInterval::freeSubIntervals()
{
    if (subIntervals)
    {
        for (auto interval : *subIntervals) {
            interval.freeSubIntervals();
        }
        delete subIntervals;
    }

}

QualQuantizer::QualQuantizer(long *nObservationsPerQual, int nLevels, int minInterestingQual)
{
    this->nObservationsPerQual = nObservationsPerQual;
    this->nLevels = nLevels;
    this->minInterestingQual = minInterestingQual;

    if (nLevels < 0)
        throw "nLevels must be >= 0";
    if (minInterestingQual < 0)
        throw "minInterestingQual must be >= 0";

    //actually run the quantizer
    this->quantizedIntervals = quantize();

    //store the map
    this->originalToQuantizeMap = intervalsToMap(quantizedIntervals);

    /*   //---for debugging
    for (int i=0; i<NQualsInHistogram; i++)
    {
        cout << i << "\t" <<  (int)originalToQuantizeMap[i] << endl;
    }   */
}

QualQuantizer::~QualQuantizer()
{
    for( auto interval : *quantizedIntervals)
    {
        interval.freeSubIntervals();
    }
    delete quantizedIntervals;
}

int QualQuantizer::getNQualsInHistogram()
{
    return NQualsInHistogram;
}

set<QualQuantizer::QualInterval > * QualQuantizer::quantize()
{
    set<QualInterval > * intervals = new set<QualInterval>;
    for (int qStart=0; qStart < NQualsInHistogram; qStart++)
    {
        long nObs = nObservationsPerQual[qStart];
        double errorRate = QualityUtils::qualToErrorProb((uint8_t)qStart);
        double nErrors = nObs * errorRate;
        intervals->emplace(qStart, qStart, nObs, (int)floor(nErrors), 0, qStart, this->minInterestingQual);
    }

    // greedy algorithm:
    // while ( n intervals >= nLevels)
    //   find intervals to merge with least penalty
    //   merge it
    while (intervals->size() > nLevels)
    {
        mergeLowestPenaltyIntervals(intervals);
    }

    return intervals;
}

void QualQuantizer::mergeLowestPenaltyIntervals(set<QualInterval> * intervals)
{
    set<QualInterval>::iterator itr1;
    set<QualInterval>::iterator itr2;

    itr1 = intervals->begin();
    itr2 = itr1;
    itr2 ++;    // skip one

    QualInterval * minMerge = nullptr;

    int lastMergeOrder = 0;
    int i=0;
    while (itr2 != intervals->end())
    {
        QualInterval* left = const_cast<QualInterval* >(&(*itr1));
        QualInterval* right = const_cast<QualInterval* >(&(*itr2));
        QualInterval* merged = left->merge(right);
        lastMergeOrder = max(max(lastMergeOrder, left->mergeOrder), right->mergeOrder);


        if (minMerge == nullptr)
        {
            minMerge = merged;
        }
        else if (merged->getPenalty() < minMerge->getPenalty())
        {
            delete minMerge->subIntervals;
            delete minMerge;
            minMerge = merged;
        }else{  // minMerge != NULL
            delete merged->subIntervals;
            delete merged;
        }

        itr1++;
        itr2++;
        i++;
    }


    set<QualInterval >::iterator itr;
    for (itr = minMerge->subIntervals->begin(); itr != minMerge->subIntervals->end(); itr++)
    {
        intervals->erase(*itr);
    }
    //cout << "inserting " << (*minMerge).qStart << "-" << (*minMerge).qEnd << endl;

    minMerge->mergeOrder = lastMergeOrder + 1;
    intervals->insert(*minMerge);
    delete minMerge;
    /*
    //---debuging   Output the whole interval set
    set<QualInterval >::iterator itr3;
    for (itr3 = intervals->begin(); itr3 != intervals->end(); itr3++)
    {
        cout << itr3->qStart << "-" << itr3->qEnd << "\t";
    }
    cout << endl;
    */

    //---debugging
    //set<QualInterval* >::iterator itr;
    /*
    cout << "subIntervals: " << endl;
    for (itr = minMerge->subIntervals->begin(); itr != minMerge->subIntervals->end(); itr++) {
        cout << (*itr).qStart << " : " << (*itr).qEnd << endl;
        cout << (*itr).subIntervals->size() << endl;
    }
    cout << "size of intervals: " << intervals->size() << endl;
    */

}

QualQuantizer::QualInterval* QualQuantizer::QualInterval::merge( QualInterval * toMerge)
{
    QualInterval* left = this->compareTo(toMerge) < 0? this : toMerge;
    QualInterval* right = this->compareTo(toMerge) < 0? toMerge : this;

    if (left->qEnd + 1 != right->qStart)
        throw "attempting to merge non-contiguous intervals";   //---???

    long nCombineObs = left->nObservations + right->nObservations;
    long nCombineErr = left->nErrors + right->nErrors;

    int level = max(left->level, right->level) + 1;
    set<QualInterval > * subIntervals = new set<QualInterval >;
    subIntervals->insert(*left);
    subIntervals->insert(*right);
    QualInterval* merged = new QualInterval(left->qStart, right->qEnd, nCombineObs, nCombineErr, level, -1, this->minQual ,subIntervals);

    return merged;
}

char * QualQuantizer::intervalsToMap(set <QualInterval> *intervals)
{
    char * map = new char[NQualsInHistogram];
    memset(map, -1, NQualsInHistogram);
    set<QualInterval>::iterator itr;
    for (itr = intervals->begin(); itr != intervals->end(); itr ++)
    {
        for (int q=itr->qStart; q<=itr->qEnd; q++)
        {
            map[q] = itr->getQual();
        }
    }

    return map;
}

char * QualQuantizer::getOriginalToQuantizedMap()
{
    return originalToQuantizeMap;
}
