cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})
project(balhbt CXX C)
message("CXX is ${CXX}")

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR}")

file(GLOB balhbt_src src/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(balhbt STATIC ${balhbt_src})
target_link_libraries(balhbt ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})
message("Libraries: ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/hp/alice/software/b3d2/include)


