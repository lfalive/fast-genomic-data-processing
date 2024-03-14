include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")
include_directories("${PROJECT_SOURCE_DIR}/src/param")

aux_source_directory(. DIR_LIB_SRCS)

add_library(HaplotypeFunctions ${DIR_LIB_SRCS} EventMap.cpp EventMap.h)