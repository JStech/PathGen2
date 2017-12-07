project(PathGen)
cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/")

string(TOLOWER ${PROJECT_NAME} LIBRARY_NAME)

# Mac ONLY. Define Relative Path on Mac OS
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

include(install_package)
include(SetPlatformVars)

set(PathGen_VERSION_MAJOR 0)
set(PathGen_VERSION_MINOR 1)
set(PACKAGE_VERSION ${PathGen_VERSION_MAJOR}.${PathGen_VERSION_MINOR})
set(PathGen_VERSION ${PACKAGE_VERSION})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()


if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Release )
endif()

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()


###############################################################################
# Set up options
option(BUILD_APPLICATIONS "Build Applications" ON)
SET(N_CORES 8 CACHE STRING "Using N number of cores for parallel build")


#####################################################################
## Find required packages
find_package(Eigen3 REQUIRED)
find_package(Glog REQUIRED)
find_package(GFlags REQUIRED)
find_package(Sophus REQUIRED)

set( REQUIRED_INCLUDE_DIRS
    ${GLOG_INCLUDE_DIRS}
    ${GFLAGS_INCLUDE_DIRS}
    ${EIGEN_INCLUDE_DIRS}
    ${Sophus_INCLUDE_DIRS}
    )

set( REQUIRED_LIBRARIES
    ${GLOG_LIBRARIES}
    ${GFLAGS_LIBRARIES}
    ${Sophus_LIBRARIES}
    )



#####################################################################
## Sources
set( INC_DIR include/${LIBRARY_NAME} )

set( HEADERS
     ${INC_DIR}/EulerSpline.h
     ${INC_DIR}/PoseSpline.h
     ${INC_DIR}/QuaternionAddons.h
     ${INC_DIR}/QuaternionSpline.h
     ${INC_DIR}/Types.h
     ${INC_DIR}/Utils.h
     ${INC_DIR}/Measurements.h
     ${INC_DIR}/IMUGenerator.h
     ${INC_DIR}/PathGenerator.h
   )

foreach(INC_FILE ${HEADERS})
  string(REPLACE "include/" "" INC_FILE ${INC_FILE})
  set(PATHGEN_INCLUDES "${PATHGEN_INCLUDES}\n#include <${INC_FILE}>")
endforeach(INC_FILE)

configure_file(
  include/pathgen/PathGen.h.in
  include/pathgen/PathGen.h
)

set( SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )

set( SOURCES
    ${SRC_DIR}/QuaternionSpline.cpp
    ${SRC_DIR}/PathGenerator.cpp
    ${SRC_DIR}/Utils.cpp
    ${SRC_DIR}/IMUGenerator.cpp
    )

list(APPEND REQUIRED_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)
list(APPEND REQUIRED_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include)

include_directories(${REQUIRED_INCLUDE_DIRS})

add_library(${LIBRARY_NAME} SHARED ${SOURCES} ${HEADERS})
set_target_properties( ${LIBRARY_NAME} PROPERTIES VERSION "${PACKAGE_VERSION}" SOVERSION "${PACKAGE_VERSION}" )
target_link_libraries(${LIBRARY_NAME} ${REQUIRED_LIBRARIES})

install_package(
    PKG_NAME ${PROJECT_NAME}
    LIB_NAME ${LIBRARY_NAME}
    VERSION ${PACKAGE_VERSION}
    DESCRIPTION "Path and IMU data generation library"
    INSTALL_INCLUDE_DIR true
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    INCLUDE_DIRS ${REQUIRED_INCLUDE_DIRS}
    LINK_LIBS ${REQUIRED_LIBRARIES}
    )
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)



#####################################################################
## Add Components
# Build application
if(BUILD_APPLICATIONS)
  MESSAGE(STATUS "BUILDING APPLICATIONS")
  set( PathGen_DIR ${CMAKE_CURRENT_BINARY_DIR} )
  add_subdirectory(applications)
endif()
