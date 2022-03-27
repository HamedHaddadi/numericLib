# determine version
cmake_minimum_required(VERSION 3.23.0-rc2)
enable_language(CXX)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()
message("Generated in ${CMAKE_BUILD_TYPE} Mode")
message(${CMAKE_CXX_COMPILER_ID})

if (${CMAKE_CXX_COMPILER_ID} STREQUAL AppleClang)
    set(CMAKE_CXX_FLAGS "-std=c++2a -g -O0")
endif()

# define project (name of the executable)
project(numLib)

include_directories("arrays")
include_directories("linearSystems")
include_directories("utils")

add_subdirectory(./matrixOperations)
