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

#add_executable(serial serial.cpp)
add_definitions(-std=c++11)
#target_link_libraries(serial clustering)

include_directories(${CMAKE_SOURCE_DIR}/../lib ${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES ( "$ENV{Eigen3_INCLUDE_DIR}" )
link_directories(${CMAKE_SOURCE_DIR}/../lib/lib/)
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#set(SOURCES serial.cpp common.cpp)
#set(HEADERS common.h )
#add_executable (serial ${SOURCES} ${HEADERS})
add_executable (serial serial.cpp)
target_link_libraries (serial clustering)

#INCLUDE_DIRECTORIES ( "${Eigen3_INCLUDE_DIR}" )
#cmake_minimum_required (VERSION 3.0)
#project (myproject)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#add_executable (serial serial.cpp)
#target_link_libraries (serial Eigen3::Eigen)
