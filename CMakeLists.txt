cmake_minimum_required(VERSION 2.8)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(PolSARSignatures)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

if ( CMAKE_COMPILER_IS_GNUCXX )
    set(MY_CXX_FLAGS  "-Wall -std=c++11 -fmessage-length=0 -v -fpermissive -fopenmp -flto -funroll-loops")
    set(CMAKE_CXX_FLAGS  "${MY_CXX_FLAGS} -O3")
endif ()

find_package(OTB  REQUIRED)
include(${OTB_USE_FILE})
message(STATUS "Found OTB: ${OTB_USE_FILE}")

set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(COMMON_INCLUDE ${PROJECT_SOURCE_DIR}/include)
set(TESTING ${PROJECT_SOURCE_DIR}/Testing)
set(UTILS ${PROJECT_SOURCE_DIR}/Utilities)

include_directories("${ITK_INCLUDE_DIRS}")
include_directories (${COMMON_INCLUDE})
link_directories(${OTB_LIBRARY_PATH})

add_subdirectory(${SRC_DIR})
add_subdirectory(${TESTING})
add_subdirectory(${UTILS})
