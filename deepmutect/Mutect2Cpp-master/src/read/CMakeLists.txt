include_directories("${PROJECT_SOURCE_DIR}/src")

aux_source_directory(. DIR_LIB_SRCS)

add_library(ReadFunctions ${DIR_LIB_SRCS} ReadCache.cpp ReadCache.h Fragment.cpp Fragment.h)

target_link_libraries(ReadFunctions PRIVATE Transformers)