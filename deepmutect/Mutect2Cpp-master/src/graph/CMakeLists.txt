include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cache")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")
include_directories("${PROJECT_SOURCE_DIR}/src/haplotype")
include_directories("${PROJECT_SOURCE_DIR}/src/path")
include_directories("${PROJECT_SOURCE_DIR}/src/read")
include_directories("${PROJECT_SOURCE_DIR}/src/smithwaterman")
aux_source_directory(. DIR_LIB_SRCS)
aux_source_directory(./set DIR_LIB_SRCS_SET)
aux_source_directory(./BaseGraph DIR_LIB_SRCS_BASEGRAPH)
aux_source_directory(./utils DIR_LIB_SRCS_UTILS)
add_library(GraphFunctions ${DIR_LIB_SRCS} ${DIR_LIB_SRCS_SET} ${DIR_LIB_SRCS_BASEGRAPH} ${DIR_LIB_SRCS_UTILS} KBestHaplotypeFinder.cpp KBestHaplotypeFinder.h KBestHaplotype.cpp KBestHaplotype.h)

# Add this flag to turn on the sort mode for comparison with GATK
#add_definitions(-DSORT_MODE)

target_link_libraries(GraphFunctions SWFunctions)