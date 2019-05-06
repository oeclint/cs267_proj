cmake_minimum_required (VERSION 3.0)
#cmake_minimum_required (VERSION 2.6) 
project (clustering)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#find_package(Eigen3 REQUIRED)
#include_directories(${Eigen3_INCLUDE_DIR})

#include_directories(${CMAKE_SOURCE_DIR}/../lib)

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
#      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
#      FORCE)
#endif(NOT CMAKE_BUILD_TYPE)

#link_directories(${CMAKE_SOURCE_DIR}/../build/lib/)

#add_executable(parallel parallel.cpp)
add_definitions(-std=c++11)
#target_link_libraries(parallel clustering)

include_directories(${CMAKE_SOURCE_DIR}/../lib ${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES ( "$ENV{Eigen3_INCLUDE_DIR}" )
link_directories(${CMAKE_SOURCE_DIR}/../lib/lib/)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#set(SOURCES parallel.cpp common.cpp)
#set(HEADERS common.h)
#add_executable (parallel ${SOURCES} ${HEADERS})
add_executable (parallel parallel.cpp)
target_link_libraries (parallel clustering)

#INCLUDE_DIRECTORIES ( "${Eigen3_INCLUDE_DIR}" )
#cmake_minimum_required (VERSION 3.0)
#project (myproject)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#add_executable (parallel parallel.cpp)
#target_link_libraries (parallel Eigen3::Eigen)
