project(balance_hbt)
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

add_custom_target(extern_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)
add_custom_target(extern_coral COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/coral/software)
add_custom_target(extern_b3d2 COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/hp/alice/software/b3d2)
add_custom_target(extern_balhbt COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/balance_hbt/software)

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/hp/alice/software/b3d2/include ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/coral/software/include ${GITHOME_MSU}/balance_hbt/software/include)        

link_directories(${GITHOME_MSU}/balance_hbt/software/lib ${GITHOME_MSU}/hp/alice/software/b3d2/lib ${GITHOME_MSU}/commonutils/software/lib ${GITHOME_MSU}/coral/software/lib)

find_library(${GITHOME_MSU}/commonutils/software/lib commonutils)
find_library(${GITHOME_MSU}/coral/software/lib coral)
find_library(${GITHOME_MSU}/hp/alice/software/b3d2/lib b3d_alice)
find_library(${GITHOME_MSU}/balance_hbt/software/lib balhbt)

add_dependencies(extern_balhbt extern_b3d2 extern_coral extern_commonutils)

#################################################
set(SOURCE_FILES ./balhbtmain.cc)
add_executable(balhbt ${CMAKE_CURRENT_SOURCE_DIR}/balhbtmain.cc)
add_dependencies(balhbt extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(balhbt ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./bwspectra.cc)
add_executable(bwspectra ${CMAKE_CURRENT_SOURCE_DIR}/bwspectra.cc)
add_dependencies(bwspectra extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(bwspectra ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

#################################################
set(SOURCE_FILES ./bw_R_CB.cc)
add_executable(bw_R_CB ${CMAKE_CURRENT_SOURCE_DIR}/bw_R_CB.cc)
add_dependencies(bw_R_CB extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(bw_R_CB ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

##################################################################################################
set(SOURCE_FILES ./protonsonly.cc)
add_executable(protonsonly ${CMAKE_CURRENT_SOURCE_DIR}/protonsonly.cc)
add_dependencies(protonsonly extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(protonsonly ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

##################################################################################################
set(SOURCE_FILES ./decaytest.cc)
add_executable(decaytest ${CMAKE_CURRENT_SOURCE_DIR}/decaytest.cc)
add_dependencies(decaytest extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(decaytest ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})

##################################################################################################
set(SOURCE_FILES ./spectratester.cc)
add_executable(spectratester ${CMAKE_CURRENT_SOURCE_DIR}/spectratester.cc)
add_dependencies(spectratester extern_commonutils)
target_link_libraries(spectratester ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ./chisquare.cc)
add_executable(chisquare ${CMAKE_CURRENT_SOURCE_DIR}/chisquare.cc)
add_dependencies(chisquare extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(chisquare ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
set(SOURCE_FILES ./balmain.cc)
add_executable(bal ${CMAKE_CURRENT_SOURCE_DIR}/balmain.cc)
add_dependencies(bal extern_balhbt extern_b3d2 extern_coral extern_commonutils)
target_link_libraries(bal ${GITHOME_MSU}/balance_hbt/software/lib/libbalhbt.a ${GITHOME_MSU}/hp/alice/software/b3d2/lib/libb3d_alice.a ${GITHOME_MSU}/coral/software/lib/libcoral.a ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
#################################################
#set(SOURCE_FILES ./crap.cc)
#add_executable(crap ${CMAKE_CURRENT_SOURCE_DIR}/crap.cc)
#add_dependencies(crap extern_commonutils)
#target_link_libraries(crap ${GITHOME_MSU}/commonutils/software/lib/libcommonutils.a ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
