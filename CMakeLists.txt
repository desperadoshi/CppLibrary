cmake_minimum_required(VERSION 3.0)

project(CppLibrary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "-std=c++11 -D USE_EIGEN3")
set(LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(HEADER_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/include")
set(BINARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

add_subdirectory("src/core")
add_subdirectory("src/app")
