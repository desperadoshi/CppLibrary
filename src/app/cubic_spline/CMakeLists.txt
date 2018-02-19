project(CubicSplineTest)

include_directories("${CMAKE_SOURCE_DIR}/src/core/cubic_spline")
link_directories("${CMAKE_SOURCE_DIR}/lib")

# Or set extra libs
# set(EXTRA_LIBS ${EXTRA_LIBS} libbrent.a)

add_executable(${PROJECT_NAME} test.cpp)

target_link_libraries(${PROJECT_NAME} cubic_spline)

add_dependencies(${PROJECT_NAME} cubic_spline)

#
# Output
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BINARY_OUTPUT_DIRECTORY}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LIBRARY_OUTPUT_DIRECTORY}")
#
# Install
#
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION "${BINARY_OUTPUT_DIRECTORY}"
  LIBRARY DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}"
  ARCHIVE DESTINATION "${LIBRARY_OUTPUT_DIRECTORY}")
# install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
#   DESTINATION "${HEADER_OUTPUT_DIRECTORY}"
#   FILES_MATCHING PATTERN "*.h")