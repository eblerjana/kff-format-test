cmake_minimum_required(VERSION 3.15)

# set the project name
project(kff-test)
project(kff-test VERSION 0.1)

set(CMAKE_BUILD_TYPE Debug)

# set binary output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# declare libraries
add_subdirectory(kff-format-api/)
add_subdirectory(src/)
