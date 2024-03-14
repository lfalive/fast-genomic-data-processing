include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/engine")
include_directories("${PROJECT_SOURCE_DIR}/src/tools")
include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")

aux_source_directory(. DIR_LIB_SRCS)
add_library(Mutect2FilteringEngineFunctions ${DIR_LIB_SRCS} ReadPositionFilter.cpp ReadPositionFilter.h MinAlleleFractionFilter.cpp MinAlleleFractionFilter.h ClusteredEventsFilter.cpp ClusteredEventsFilter.h GermlineFilter.cpp GermlineFilter.h MultiallelicFilter.cpp MultiallelicFilter.h FragmentLengthFilter.cpp FragmentLengthFilter.h PolymeraseSlippageFilter.cpp PolymeraseSlippageFilter.h)

target_link_libraries(Mutect2FilteringEngineFunctions ToolsFunctions)

