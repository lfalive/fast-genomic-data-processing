//
// Created by lhh on 5/30/22.
//

#ifndef MUTECT2CPP_MASTER_SUBSETTEDLIKELIHOODMATRIX_H
#define MUTECT2CPP_MASTER_SUBSETTEDLIKELIHOODMATRIX_H

/**
 * Fast wrapper for a LikelihoodMatrix that uses only a subset of alleles.  Useful for model comparison of different
 * allele subsets without having to copy the underlying likelihoods.
 */

#include "utils/genotyper/AlleleLikelihoods.h"
#include "Mutect2Utils.h"

template<typename EVIDENCE, typename A>
class SubsettedLikelihoodMatrix {
private:
    SampleMatrix<EVIDENCE, A>* matrix;
    shared_ptr<vector<shared_ptr<Allele>>> alleles;
    map<int, int> newToOldIndexMap;

public:
    SubsettedLikelihoodMatrix(SampleMatrix<EVIDENCE,A>* matrix, shared_ptr<vector<shared_ptr<Allele>>> alleles):matrix(matrix), alleles(alleles)
    {
        int IndexSize = alleles->size();
        for(int i=0; i<IndexSize; i++)
        {
            int newIndex =  matrix->indexOfAllele(alleles->operator[](i));
            assert(newIndex != -1);
            newToOldIndexMap.insert({i, newIndex});
        }

    }

    static std::shared_ptr<SubsettedLikelihoodMatrix<EVIDENCE,A>> excludingAllele(SampleMatrix<EVIDENCE,A>* matrix, shared_ptr<Allele> excludedAllele)
    {
        auto alleles = make_shared<vector<shared_ptr<A>>>();
        for(auto& a : matrix->alleles())
        {
            if(!basesMatch(a, excludedAllele))
                alleles->emplace_back(a);
        }
        Mutect2Utils::validateArg(alleles->size() == matrix->numberOfAlleles()-1, "More than one allele excluded.");
        return make_shared<SubsettedLikelihoodMatrix<EVIDENCE,A>>(matrix, alleles);

    }

    static bool basesMatch(shared_ptr<Allele> a, shared_ptr<Allele> b)
    {
        assert(a != nullptr && b != nullptr);
        if(a.get() == b.get())
            return true;
        if(a->getBases().get() == b->getBases().get())
            return true;

        if(a->getBasesLength() != b->getBasesLength())
            return false;
        else{
            auto bases_ = a->getBases();
            auto other_ = b->getBases();
            for(int i = 0; i < a->getBasesLength(); i++) {
                if(bases_[i] != other_[i])
                    return false;
            }
            return true;
        }
    }

    int numberOfAlleles() { return alleles->size(); }

    int evidenceCount() {
        return matrix->evidenceCount();
    }

    double get(int alleleIndex, int evidenceIndex)
    {
        assert(newToOldIndexMap.find(alleleIndex) != newToOldIndexMap.end());
        return matrix->get(newToOldIndexMap[alleleIndex], evidenceIndex);
    }

    shared_ptr<Allele> getAllele(int alleleIndex){
        return alleles->operator[](alleleIndex);
    }

    shared_ptr<vector<shared_ptr<Allele>>> getAlleles()
    {
        return alleles;
    }
};


#endif //MUTECT2CPP_MASTER_SUBSETTEDLIKELIHOODMATRIX_H
