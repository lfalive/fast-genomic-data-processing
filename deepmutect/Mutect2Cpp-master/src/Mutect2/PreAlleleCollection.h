//
// Created by lhh on 5/27/22.
//

#ifndef MUTECT2CPP_MASTER_PREALLELECOLLECTION_H
#define MUTECT2CPP_MASTER_PREALLELECOLLECTION_H

#include "parallel_hashmap/phmap.h"
#include <memory>
#include <cassert>
#include "Allele.h"
#include "Mutect2Utils.h"


template <typename X>
class PreAlleleCollection {
public:
    enum Type{
        ALT_ONLY, REF_AND_ALT
    };

private:
    std::shared_ptr<Allele> refAllele;
    X refValue;
    phmap::flat_hash_map<Allele*, X> altAlleleValueMap;   // map of Allele* tp value, you should make sure Allele* is a valid pointer
    Type type;

public:



    explicit PreAlleleCollection(Type type) : refAllele(nullptr), refValue(0), type(type)
    {

    }

    /**
     * Take an allele, REF or ALT, and update its value appropriately
     *
     * @param allele : REF or ALT allele
     * @param value :
     */
     void set(std::shared_ptr<Allele> allele, X value)
     {
         assert(allele != nullptr);

         assert(type == REF_AND_ALT || allele->getIsNonReference());
         if(allele->getIsReference())
             setRef(allele, value);
         else
             setAlt(allele, value);
     }

     void setRef(std::shared_ptr<Allele> allele, X value)
     {
         assert(allele != nullptr);
         assert(allele->getIsReference());
         Mutect2Utils::validateArg(refAllele == nullptr, "Resetting the reference allele not permitted");
        refAllele = allele;
        refValue = value;
     }

     void setAlt(std::shared_ptr<Allele> allele, X value)
     {
         assert(allele != nullptr);
         Mutect2Utils::validateArg(allele->getIsNonReference(), "Setting reference allele as alt");
         altAlleleValueMap.template insert({allele.get(), value});
     }

     /**
       * Get the value for an allele, REF or ALT
       * @param allele
       */
     X get(std::shared_ptr<Allele> allele){
         assert(allele != nullptr);
         if(allele->getIsReference())
         {
             Mutect2Utils::validateArg(*allele == *refAllele, "Requested ref allele does not match the stored ref allele");
             return getRef();
         } else {
             return getAlt(allele);
         }
     }

     X getRef()
     {
         if(type == ALT_ONLY)
             throw "Collection does not hold the REF allele";

         if(refAllele != nullptr)
         {
             return refValue;
         } else {
             throw "Collection's ref allele has not been set yet";
         }
     }

     X getAlt(std::shared_ptr<Allele> allele)
     {
         assert(allele != nullptr);
         Mutect2Utils::validateArg(allele->getIsNonReference(), "allele is not an alt allele");
         Mutect2Utils::validateArg(altAlleleValueMap.find(allele.get()) != altAlleleValueMap.end(), "Requested alt allele is not in the collection");
         return altAlleleValueMap.at(allele.get());
     }

     std::shared_ptr<std::vector<double>> asDoubleArray(std::shared_ptr<std::vector<std::shared_ptr<Allele>>> allelesInOrder)
     {
         return asDoubleArray(*allelesInOrder);
     }

     std::shared_ptr<std::vector<double>> asDoubleArray(std::vector<std::shared_ptr<Allele>>& allelesInOrder)
     {
         auto result = std::make_shared<std::vector<double>>();
         for(auto& a : allelesInOrder)
         {
             result->template emplace_back((double)get(a));
         }
         return result;
     }
};


#endif //MUTECT2CPP_MASTER_PREALLELECOLLECTION_H
