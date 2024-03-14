//
// Created by 梦想家xixi on 2021/12/6.
//

#ifndef MUTECT2CPP_MASTER_VARIANTCONTEXTBUILDER_H
#define MUTECT2CPP_MASTER_VARIANTCONTEXTBUILDER_H

#include <string>
#include <vector>
#include <set>
#include "parallel_hashmap/phmap.h"
#include "Allele.h"
#include "variantcontext/GenoTypesContext.h"
#include "variantcontext/AttributeValue.h"

class VariantContextBuilder {
private:
    // required fields
    bool fullyDecoded = false;
    std::string source;
    std::string contig;
    long start = -1;
    long stop = -1;
    std::shared_ptr<std::vector<std::shared_ptr<Allele>>> alleles;

    // optional -> these are set to the appropriate default value
    std::string ID = ".";
    std::shared_ptr<GenoTypesContext> genotypes;
    double log10PError;
    std::set<std::string> * filters;
    std::shared_ptr<std::map<std::string, AttributeValue>> attribute;
    bool attributesCanBeModified = false;
    bool filtersCanBeModified = false;
    std::set<Validation> toValidate;

public:
    VariantContextBuilder(std::string & source, std::string & contig, long start, long stop, std::shared_ptr<std::vector<std::shared_ptr<Allele>>>   alleles);

    explicit VariantContextBuilder(std::shared_ptr<VariantContext> & parent);

    /**
     * Create an empty VariantContextBuilder where all values adopt their default values.  Note that
     * source, chr, start, stop, and alleles must eventually be filled in, or the resulting VariantContext
     * will throw an error.
     */
    VariantContextBuilder();

    ~VariantContextBuilder();



    std::shared_ptr<VariantContext> make(bool leaveModifyableAsIs);

    std::shared_ptr<VariantContext> make();

    void setStart(long start);

    void setStop(long stop);

    VariantContextBuilder* setAlleles(const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> & alleles);

    /**
    * Tells this builder to use this collection of alleles for the resulting VariantContext
    *
    * @param alleles an hash set of alleles to set as the alleles of this builder
    */
    VariantContextBuilder* setAlleles(phmap::flat_hash_set<std::shared_ptr<Allele>>& alleles);

    /**
     * Tells us that the resulting VariantContext should have source field set to source
     * @param source string describing the source of the variant
     * @return this builder
    */     //---remember to free it
    VariantContextBuilder* setSource(std::string& source);

    VariantContextBuilder* setId(std::string& ID);

    VariantContextBuilder* setLoc(std::string& contig, long start, long stop);

    VariantContextBuilder* setGenotypes(std::shared_ptr<GenoTypesContext> genotypesContext);

    VariantContextBuilder* setGenotypes(std::vector<std::shared_ptr<Genotype>>& genotypes);

    VariantContextBuilder* setLog10PError(double log10PError);

    /**
   * This builder's filters are set to this value
   *
   * filters can be <code>null</code> -&gt; meaning there are no filters
   *
   * @param filters Set of strings to set as the filters for this builder
   *                This set will be copied so that external set can be
   *                safely changed.
   * @return this builder
   */
    VariantContextBuilder* setFilters(std::set<std::string> * filters);

    VariantContextBuilder* unfiltered();

    VariantContextBuilder* setAttributes(std::shared_ptr<std::map<std::string, AttributeValue>> attribute);

    VariantContextBuilder* setAttributes(std::shared_ptr<std::map<std::string, std::vector<double>>> attribute);

    VariantContextBuilder* setAttribute(const std::string& key, int value);

    VariantContextBuilder* setAttribute(const std::string& key, std::vector<double> value);

    VariantContextBuilder* setAttribute(const std::string& key, std::shared_ptr<std::vector<double>> value);
};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXTBUILDER_H
