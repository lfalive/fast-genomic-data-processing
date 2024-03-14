//
// Created by 梦想家xixi on 2021/11/1.
//

#include <cassert>
#include "MathUtils.h"
#include <cmath>
#include <cfloat>
using namespace std;


double MathUtils::digamma(int n) {
    return DIGAMMA_CACHE().get(n);
}

double MathUtils::log10ToLog(double log10) {
    return log10 * std::log(10);
}

double MathUtils::log10Factorial(int n) {
    return LOG_10_FACTORIAL_CACHE().get(n);
}

double MathUtils::fastBernoulliEntropy(const double p) {
    double product = p * (1 - p);
    return product * (11 + 33 * product) / (2 + 20 * product);
}

double MathUtils::normalDistribution(double mean, double sd, double x)
{
    assert(sd >= 0);
    assert((!isnan(mean)) && (!isinf(mean)) && (!isnan(sd)) && (!isinf(sd)) && !isnan(x) && (!isinf(x)));
    return exp(-(x - mean) * (x - mean) / (2.0 * sd * sd)) / (sd * sqrt(2.0 * M_PI));
}

std::vector<double> * MathUtils::normalizeSumToZero( std::vector<double> * array)
{
    assert(array != nullptr);
    if(array->empty())
        return array;

    double sum = 0.0;
    for(auto &element: *array)
    {
        sum += element;
    }
    for(auto &element: *array)
    {
        element = element/sum;
    }
    return array;
}

DigammaCache &MathUtils::DIGAMMA_CACHE() {
    static thread_local DigammaCache digamma;
    return digamma;
}

Log10FactorialCache &MathUtils::LOG_10_FACTORIAL_CACHE() {
    static thread_local Log10FactorialCache log10Factorial;
    return log10Factorial;
}


int MathUtils::fastRound(double d)
{
    return d > 0.0 ? (int)(d+0.5) : (int)(d-0.5);
}

double MathUtils::log10BinomialProbability(int n, int k, double log10p)
{
    if (log10p == std::numeric_limits<double>::infinity())
    {
        return k == 0 ? 0 : std::numeric_limits<double>::infinity();
    }
    double log10OneMinusP = log10(1 - pow(10.0, log10p));
    return log10BinomialCoefficient(n, k) + log10p * k + log10OneMinusP * (n - k);
}

double MathUtils::log10Gamma(double x)
{
    return lgamma(x) * M_LOG10E;
}

double MathUtils::log10Fractorial(int n)
{
    return LOG_10_FACTORIAL_CACHE().get(n);
}

double MathUtils::log10BinomialCoefficient(int n, int k)
{
    assert(n >= 0);
    assert(k <= n && k >= 0);
    return log10Fractorial(n) - log10Fractorial(k) - log10Fractorial(n-k);
}

double MathUtils::sum(vector<double> &values) {
    double sum = accumulate(values.begin(), values.end(), 0.0);
    return sum;
}

double MathUtils::sum(shared_ptr<vector<double>> values) {
    double sum = accumulate(values->begin(), values->end(), 0.0);
    return sum;
}

shared_ptr<vector<double>> MathUtils::applyToArray(shared_ptr<vector<double>> array, double (*func)(double)) {
    return applyToArray(*array, func);
}

shared_ptr<vector<double>> MathUtils::applyToArray(vector<double> &array, double (*func)(double)) {
    assert(!array.empty());
    assert(func != nullptr);
    auto result = make_shared<vector<double>>(array.size());
    for(int m=0; m<result->size(); m++)
        result->operator[](m) = func(array[m]);
    return result;
}

shared_ptr<vector<double>> MathUtils::applyToArrayInPlace(shared_ptr<vector<double>> array, function<double(double)> func) {
    assert(array != nullptr);
    assert(func != nullptr);
    applyToArrayInPlace(*array, func);
    return array;
}

void MathUtils::applyToArrayInPlace(vector<double> &array, function<double(double)> func) {
    for(double & m : array)
    {
        m = func(m);
    }
}

shared_ptr<vector<double>>
MathUtils::sumArrayFunction(int min, int max, function<shared_ptr<vector<double>>(int)> func) {
    Mutect2Utils::validateArg(max >= min, "max must be at least as great as min");
    auto result = func(min);
    for(int n = min+1; n < max; n++)
    {
        auto newValues = func(n);
        Mutect2Utils::validateArg(newValues->size() == result->size(), "array function returns different sizes for different inputs!");
        for(int i=0; i<result->size(); i++)
        {
            result->operator[](i) += newValues->operator[](i);
        }
    }
    return result;
}

int MathUtils::maxElementIndex(shared_ptr<vector<double>> array) {
    assert(array != nullptr);
    return maxElementIndex(*array);
}

int MathUtils::maxElementIndex(vector<double>& array) {
    int size = array.size();
    int maxI = 0;
    for(int i=1; i<size; i++)
    {
        if(array[i] > array[maxI])
            maxI = i;
    }
    return maxI;
}

double MathUtils::distance1(vector<double> &p1, vector<double> &p2) {
    double sum = 0.0;
    for(int i=0; i<p1.size(); i++)
    {
        sum += abs(p1[i] - p2[i]);
    }
    return sum;
}

shared_ptr<vector<double>> MathUtils::normalizeSumToOne(shared_ptr<vector<double>> array) {
    assert(array != nullptr);
    if(array->empty())
        return array;

    double sum = MathUtils::sum(array);
    Mutect2Utils::validateArg(sum >= 0.0, "Values in probability array sum to a negative number ");
    return applyToArrayInPlace(array, [sum](double x){return x / sum;});
}

int MathUtils::median(std::vector<int> &values) {
    assert(!values.empty());
    sort(values.begin(), values.end());
    int size = values.size();
    if(size % 2 != 0)
        return values[size/2];
    double median = (values[(size-1)/2] + values[size/2]) / 2.0;
    return round(median);
}

std::vector<double> MathUtils::normalizeLog10(vector<double> &array, bool normalizeLog10, bool inPlace) {
    double log10Sum = log10SumLog10(array);
    std::vector<double> res = std::vector<double>(array.size());
    if(inPlace) {
        for(int i = 0; i < array.size(); i++) {
            array[i] = array[i] - log10Sum;
            res[i] = array[i];
        }
    } else {
        for(int i = 0; i < array.size(); i++) {
            res[i] = array[i] - log10Sum;
        }
    }
    if(normalizeLog10) {
        return res;
    } else {
        applyToArrayInPlace(res, [](double x){
            return std::pow(10.0, x);
        });
        return res;
    }
}

double MathUtils::log10SumLog10(const std::vector<double>& log10Values, int start, int finish) {
    if (start >= finish) {
        return -DBL_MAX;
    }
    int maxIndex = maxElementIndex(log10Values, start, finish);
    double maxValue = log10Values[maxIndex];
    if(maxValue == -DBL_MAX) {
        return maxValue;
    }
    double sum = 1.0;
    for (int i = start; i < finish; i++) {
        double curVal = log10Values[i];
        if (i == maxIndex || curVal == -DBL_MAX) {
            continue;
        } else {
            double scaled_val = curVal - maxValue;
            sum += std::pow(10.0, scaled_val);
        }
    }
    if ( isinf(sum) || sum == DBL_MAX ) {
        throw std::invalid_argument("log10 p: Values must be non-infinite and non-NAN");
    }
    return maxValue + (sum != 1.0 ? std::log10(sum) : 0.0);
}

int MathUtils::maxElementIndex(const vector<double> &array, int start, int endIndex) {
    int maxI = start;
    for (int i = (start+1); i < endIndex; i++) {
        if (array[i] > array[maxI])
            maxI = i;
    }
    return maxI;
}

double MathUtils::log10SumLog10(const vector<double> &log10Values) {
    return log10SumLog10(log10Values, 0, log10Values.size());
}

double MathUtils::logToLog10(const double ln) {
    return ln * 0.4342944819032518;
}
