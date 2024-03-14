include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")

aux_source_directory(. DIR_LIB_SRCS)

add_library(SWFunctions ${DIR_LIB_SRCS})

target_link_libraries(SWFunctions CigarFunctions)