cmake_minimum_required(VERSION 3.15)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# set the project name
project(kff-test)
project(kff-test VERSION 0.1)

set(CMAKE_BUILD_TYPE Debug)

# set binary output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# declare libraries
add_subdirectory(kff-format-api-cpp/)
add_subdirectory(src/)

