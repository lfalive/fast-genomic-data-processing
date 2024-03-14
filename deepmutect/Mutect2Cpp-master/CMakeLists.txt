cmake_minimum_required(VERSION 3.17)
project(Mutect2Cpp_master)

set(CMAKE_CXX_STANDARD 17)

#libtorch required
set(CMAKE_PREFIX_PATH "/usr/libtorch/")
find_package(Torch REQUIRED)

#compile flag
set(CMAKE_CXX_FLAGS "-g -pthread")
set(CMAKE_CXX_FLAGS_RELEASE -O3)
#set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=leak")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

#linker flag
#set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=leak")

# counterpart of gcc -I
include_directories(${PROJECT_SOURCE_DIR}/../htslib)
include_directories(src)

aux_source_directory(./src DIR_SRCS)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cache")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")
include_directories("${PROJECT_SOURCE_DIR}/src/haplotype")
include_directories("${PROJECT_SOURCE_DIR}/src/path")
include_directories("${PROJECT_SOURCE_DIR}/src/read")
include_directories("${PROJECT_SOURCE_DIR}/src/graph")
include_directories("${PROJECT_SOURCE_DIR}/src/intel/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/intel/pairhmm")
include_directories("${PROJECT_SOURCE_DIR}/src/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/Mutect2")
include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")
include_directories("${PROJECT_SOURCE_DIR}/src/annotator")
include_directories("${PROJECT_SOURCE_DIR}/src/filtering")
include_directories("${PROJECT_SOURCE_DIR}/src/clustering")
include_directories("${PROJECT_SOURCE_DIR}/src/basicshortmutpileup")

add_subdirectory(./src/cache)
add_subdirectory(./src/cigar)
add_subdirectory(./src/haplotype)
add_subdirectory(./src/path)
add_subdirectory(./src/transformers)
add_subdirectory(./src/read)
add_subdirectory(./src/graph)
add_subdirectory(./src/intel/smithwaterman)
add_subdirectory(./src/intel/pairhmm)
add_subdirectory(./src/smithwaterman)
add_subdirectory(./src/Mutect2)
add_subdirectory(./src/variantcontext)
add_subdirectory(./src/haplotypecaller)
add_subdirectory(./src/utils)
add_subdirectory(./src/tools)
add_subdirectory(./src/samtools)
add_subdirectory(./src/clipping)
add_subdirectory(./src/annotator)
add_subdirectory(src/filtering)
add_subdirectory(./src/clustering)
add_subdirectory(./src/basicshortmutpileup)

add_executable(Mutect2Cpp_master ${DIR_SRCS} src/model/model.cpp src/model/model.h src/engine/ReferenceContext.cpp src/engine/ReferenceContext.h src/engine/AlignmentContext.cpp src/engine/AlignmentContext.h src/pileRead.h
        src/trie/trieNode.cpp src/trie/trieNode.h src/trie/buildTreeUtils.cpp src/trie/buildTreeUtils.h)

# export MALLOC_CONF=background_thread:true,dirty_decay_ms:40000,metadata_thp:always,tcache_max:65536
target_link_libraries(Mutect2Cpp_master libjemalloc.a -ldl)

# Add this flag to turn on the sort mode for comparison with GATK
#add_definitions(-DSORT_MODE)

# counterpart of gcc -L
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_libraries(Mutect2Cpp_master ${PROJECT_SOURCE_DIR}/../htslib/libhts.a -lz -llzma -lbz2 -ldl)
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    target_link_libraries(Mutect2Cpp_master ${PROJECT_SOURCE_DIR}/../htslib/libhts.dylib)
endif ()
target_link_libraries(Mutect2Cpp_master CigarFunctions)
target_link_libraries(Mutect2Cpp_master ReadFunctions)
target_link_libraries(Mutect2Cpp_master CacheFunctions)
target_link_libraries(Mutect2Cpp_master GraphFunctions)
target_link_libraries(Mutect2Cpp_master HaplotypeFunctions)
target_link_libraries(Mutect2Cpp_master PathFunctions)
target_link_libraries(Mutect2Cpp_master smithwatermanFunctions)
target_link_libraries(Mutect2Cpp_master SWFunctions)
target_link_libraries(Mutect2Cpp_master Mutect2Functions)
target_link_libraries(Mutect2Cpp_master VariantContext)
target_link_libraries(Mutect2Cpp_master HaplotypeCallerFunctions)
target_link_libraries(Mutect2Cpp_master UtilsFunctions)
target_link_libraries(Mutect2Cpp_master ToolsFunctions)
target_link_libraries(Mutect2Cpp_master SamtoolsFunctions)
target_link_libraries(Mutect2Cpp_master ClippingFunctions)
target_link_libraries(Mutect2Cpp_master gkl_pairhmm)
target_link_libraries(Mutect2Cpp_master "${TORCH_LIBRARIES}")
target_link_libraries(Mutect2Cpp_master AnnotatorFunctions)
target_link_libraries(Mutect2Cpp_master Mutect2FilteringEngineFunctions)
target_link_libraries(Mutect2Cpp_master clusteringFunctions)
target_link_libraries(Mutect2Cpp_master basicshortmutpileupFunctions)


