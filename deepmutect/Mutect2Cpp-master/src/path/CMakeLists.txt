include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/graph")
include_directories("${PROJECT_SOURCE_DIR}/src/cache")

aux_source_directory(. DIR_LIB_SRCS)

add_library(PathFunctions ${DIR_LIB_SRCS})
