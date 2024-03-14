include_directories("${PROJECT_SOURCE_DIR}/src")

aux_source_directory(. DIR_LIB_SRCS)

add_library(CigarFunctions ${DIR_LIB_SRCS} TextCigarCodec.cpp TextCigarCodec.h)
