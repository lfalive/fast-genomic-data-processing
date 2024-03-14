

file(GLOB RECALIBRATION_SOURCE "./*.cpp")
add_library(Recalibration STATIC ${RECALIBRATION_SOURCE}  ${PROJECT_SOURCE_DIR}/src/MathUtils.cpp)