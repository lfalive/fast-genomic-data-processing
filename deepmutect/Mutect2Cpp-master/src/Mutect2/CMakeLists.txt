include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cache")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")
include_directories("${PROJECT_SOURCE_DIR}/src/haplotype")
include_directories("${PROJECT_SOURCE_DIR}/src/path")
include_directories("${PROJECT_SOURCE_DIR}/src/read")
include_directories("${PROJECT_SOURCE_DIR}/src/graph")
include_directories("${PROJECT_SOURCE_DIR}/src/intel/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/param")

aux_source_directory(. DIR_LIB_SRCS)

add_library(Mutect2Functions ${DIR_LIB_SRCS} AssemblyResultSet.cpp AssemblyResultSet.h SomaticGenotypeEngine.cpp SomaticGenotypeEngine.h PreAlleleCollection.cpp PreAlleleCollection.h SomaticLikelihoodsEngine.cpp SomaticLikelihoodsEngine.h SubsettedLikelihoodMatrix.h)

# Add this flag to turn on the sort mode for comparison with GATK
#add_definitions(-DSORT_MODE)

target_link_libraries(Mutect2Functions HaplotypeFunctions)
target_link_libraries(Mutect2Functions PathFunctions)
target_link_libraries(Mutect2Functions GraphFunctions)
target_link_libraries(Mutect2Functions ReadFunctions)