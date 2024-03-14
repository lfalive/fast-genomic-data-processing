include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")
include_directories("${PROJECT_SOURCE_DIR}/src")

aux_source_directory(. DIR_LIB_SRCS)

add_library(VariantContext ${DIR_LIB_SRCS} VCFConstants.h GenoTypesContext.cpp GenoTypesContext.h CommonInfo.cpp CommonInfo.h GenotypeLikelihoods.cpp GenotypeLikelihoods.h GenotypeLikelihoodsAllelePair.cpp GenotypeLikelihoodsAllelePair.h utils/GeneralUtils.cpp utils/GeneralUtils.h GenotypeType.h FastGenotype.h builder/GenotypeBuilder.cpp builder/GenotypeBuilder.h builder/VariantContextBuilder.cpp builder/VariantContextBuilder.h VariantContextUtils.cpp VariantContextUtils.h AttributeValue.cpp AttributeValue.h VCFWriter.h VCFWriter.cpp)