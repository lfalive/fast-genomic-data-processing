include_directories("${PROJECT_SOURCE_DIR}/src")

aux_source_directory(. DIR_LIB_SRCS)
add_library(basicshortmutpileupFunctions ${DIR_LIB_SRCS} CombinatoricsUtils.cpp CombinatoricsUtils.h ArithmeticUtils.cpp ArithmeticUtils.h)