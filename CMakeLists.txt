cmake_minimum_required(VERSION 3.0)

project(CubicSpline)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

### Find Eigen3
set(EIGEN3_ROOT "${CMAKE_SOURCE_DIR}/contrib/eigen3/include/eigen3")
find_package (Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
# set(LIBS ${LIBS} ${Eigen3_LIBRARIES})

include_directories("src/core")

file(GLOB_RECURSE ALL_SRC
    "src/*.cpp"
)

set(CMAKE_CXX_FLAGS "-std=c++11")

add_executable(${PROJECT_NAME} ${ALL_SRC})
