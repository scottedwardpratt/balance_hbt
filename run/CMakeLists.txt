project(balance_hbt CXX)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)
add_custom_target(extern_b3d2 COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/hp/alice/software/b3d2)
add_custom_target(extern_balhbt COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/balance_hbt/software)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3 -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-result -std=c++11 -mfpmath=sse")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/hp/alice/software/b3d2/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/balance_hbt/software/include)        

link_directories(${GITHOME_MSU}/balance_hbt/software/lib ${GITHOME_MSU}/hp/alice/software/b3d2/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)
find_library(${GITHOME_MSU}/coral/software/lib coral)
find_library(${GITHOME_MSU}/hp/alice/software/b3d2/lib b3d_alice)
find_library(${GITHOME_MSU}/balance_hbt/software/lib balhbt)

#################################################

set(SOURCE_FILES ./balance_hbt.cc)
SET_SOURCE_FILES_PROPERTIES(
      balance_hbt.cc PROPERTIES OBJECT_DEPENDS
    ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a
 )
add_executable(balance_hbt ${SOURCE_FILES})
add_dependencies(balance_hbt extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(balance_hbt balhbt b3d_alice coral commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
