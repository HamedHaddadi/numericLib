cmake_minimum_required(VERSION 3.23.0-rc2)
project(sanity)
enable_language(CXX)

set (EXECUTABLE_NAME "sanityCheck")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

message("Generated in ${CMAKE_BUILD_TYPE} Mode")
message(${CMAKE_CXX_COMPILER_ID})

if (${CMAKE_CXX_COMPILER_ID} STREQUAL AppleClang)
    set(CMAKE_CXX_FLAGS "-std=c++2a -g -O0")
endif()

include_directories(../arrays)
include_directories(../linearSystems)
include_directories(../utils)

# define project (name of the executable)
add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")