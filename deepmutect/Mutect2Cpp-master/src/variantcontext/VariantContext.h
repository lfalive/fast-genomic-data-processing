//
// Created by lhh on 11/11/21.
//

#ifndef MUTECT2CPP_MASTER_VARIANTCONTEXT_H
#define MUTECT2CPP_MASTER_VARIANTCONTEXT_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "htslib/sam.h"
#include "haplotype/Allele.h"
#include "VCFConstants.h"
#include "GenoTypesContext.h"
#include "CommonInfo.h"

class GenoTypesContext;

enum VariantContextType {
	VariantContext_NO_VARIATION,
	VariantContext_SNP,
	VariantContext_MNP,
	VariantContext_INDEL,
	VariantContext_SYMBOLIC,
	VariantContext_MIXED,
	VariantContext_NULL
};

enum Validation {
	ALLELES,
	GENOTYPES,
	FILTERS
};

class VariantContext {
private:
	std::string contig;
	hts_pos_t start;
	hts_pos_t stop;
	std::string ID;
	std::shared_ptr<Allele> REF;
	std::shared_ptr<Allele> ALT;
    std::vector<int> filters;
	CommonInfo commonInfo;
	bool fullyDecoded;

	/** A set of the alleles segregating in this context */

	void validateStop();

	bool validate(const std::set<Validation> &validationToPerform);

	static std::vector<std::shared_ptr<Allele>> makeAlleles(std::vector<std::shared_ptr<Allele>> &alleles);

	void validateAlleles();

	void validateGenotypes();

	void determineType();

	void determinePolymorphicType();

	static VariantContextType
	typeOfBiallelicVariant(const std::shared_ptr<Allele> &ref, const std::shared_ptr<Allele> &allele);

protected:
    VariantContextType type;
    std::vector<std::shared_ptr<Allele>> alleles;
    std::shared_ptr<GenoTypesContext> genotypes;


public:
    VariantContext(std::string &source,
                   std::string &ID,
                   std::string &contig,
                   long start,
                   long stop,
                   const std::shared_ptr<std::vector<std::shared_ptr<Allele>>> &alleles,
                   std::shared_ptr<GenoTypesContext> genotypes,
                   double log10PError,
                   std::set<std::string> *filters,  std::shared_ptr<std::map<std::string, AttributeValue>> attributes,
                   bool fullyDecoded,
                   const std::set<Validation> &validationToPerform
                   );

    ~VariantContext();

	/**
	 * the actual constructor.  Private access only
	 *
	 * @param source          source
	 * @param contig          the contig
	 * @param start           the start base (one based)
	 * @param stop            the stop reference base (one based)
	 * @param alleles         alleles
	 * @param genotypes       genotypes map
	 * @param log10PError  qual
	 * @param filters         filters: use null for unfiltered and empty set for passes filters
	 * @param attributes      attributes
	 * @param validationToPerform     set of validation steps to take
	 */

	bool hasAttribute(const std::string &key);

	int getAttributeAsInt(const std::string &key, int defaultValue);

	std::vector<int> getAttributeAsIntVector(const std::string &key, std::vector<int> defaultValue);

	int getEnd() const;

	int getStart();

	bool isBiallelic();

	int getNAlleles();

	bool isSNP();

	bool isSimpleDeletion();

	bool isSimpleInsertion();

	bool isSimpleIndel();

	bool isIndel();

	VariantContextType getType();

	bool hasSymbolicAlleles();

	std::vector<std::shared_ptr<Allele>> &getAlleles();

	static bool hasSymbolicAlleles(const std::vector<std::shared_ptr<Allele>> &alleles);

	std::shared_ptr<Allele> getReference();

	bool hasAllele(const std::shared_ptr<Allele> &allele);

	bool hasAllele(const std::shared_ptr<Allele> &allele, bool ignoreRefState);

	bool hasAllele(const std::shared_ptr<Allele> &allele, bool ignoreRefState, bool considerRefAllele);

	std::vector<std::shared_ptr<Allele>> getAlternateAlleles();

	std::shared_ptr<Allele> getAlternateAllele(int i);

	const std::map<std::string, AttributeValue> &getAttributes();

	std::shared_ptr<std::map<std::string, AttributeValue>> getAttributesAsPointer();

	std::string &getContig();

	int getContigInt();

	std::set<std::string>& getFilter();

	std::set<std::string> *getFiltersMaybeNull();

	std::shared_ptr<GenoTypesContext> getGenotypes();

	std::string &getID();

	bool hasID();

	double getLog10PError();

	std::string &getSource();

	bool isFullyDecoded() const;

	std::string getTypeString();

	bool isNotFiltered();

	bool isFiltered();

	/**
    * convenience method for variants
    *
    * @return true if this is a variant allele, false if it's reference
    */
	bool isVariant();

	bool filtersWereApplied();

	/**
    * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
    *
    * @return chromosome count
    */
	int getCalledChrCount();

	/**
     * Returns the number of chromosomes carrying any allele in the genotypes (i.e., excluding NO_CALLS)
     *
     * @param sampleIds IDs of samples to take into account. If empty then all samples are included.
     * @return chromosome count
     */
	int getCalledChrCount(std::set<std::string>& sampleIds);

	/**
     * @return true if the context has associated genotypes
     */
	bool hasGenotypes();

	/**
   * @return the number of samples in the context
   */
   int getNSamples();

   std::vector<int> getIndelLengths();

   void addFilter(int index);

   std::vector<int> getFilters();
};


#endif //MUTECT2CPP_MASTER_VARIANTCONTEXT_H
