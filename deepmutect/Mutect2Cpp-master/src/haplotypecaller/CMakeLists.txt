include_directories("${PROJECT_SOURCE_DIR}/src/intel/pairhmm")

aux_source_directory(. DIR_LIB_SRCS)
add_library(HaplotypeCallerFunctions ${DIR_LIB_SRCS} ReferenceConfidenceModel.cpp ReferenceConfidenceModel.h PairHMMLikelihoodCalculationEngine.cpp PairHMMLikelihoodCalculationEngine.h LikelihoodEngineArgumentCollection.cpp LikelihoodEngineArgumentCollection.h PairHMMNativeArgumentCollection.h CalledHaplotypes.cpp CalledHaplotypes.h LocationAndAlleles.cpp LocationAndAlleles.h ReadForPairHMM.cpp ReadForPairHMM.h)

# Add this flag to turn on the sort mode for comparison with GATK
#add_definitions(-DSORT_MODE)

target_link_libraries(HaplotypeCallerFunctions VariantContext)
target_link_libraries(HaplotypeCallerFunctions gkl_pairhmm)