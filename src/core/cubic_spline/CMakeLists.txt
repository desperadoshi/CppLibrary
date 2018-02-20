project(cubic_spline)

include_directories("${CMAKE_SOURCE_DIR}/src/core/misc")

### Find Eigen3
set(EIGEN3_ROOT "${CMAKE_SOURCE_DIR}/contrib/eigen3/include/eigen3")
find_package (Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories("${CMAKE_SOURCE_DIR}/src/core/linalg")

add_library(${PROJECT_NAME} STATIC cubic_spline.cpp)

#
# Output
#
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}")
#
# Install
#
install(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}"
  ARCHIVE DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  DESTINATION "${HEADER_OUTPUT_DIRECTORY}"
  FILES_MATCHING PATTERN "*.h")
