project(CubicSplineTest)

include_directories("${CMAKE_SOURCE_DIR}/src/core/cubic_spline")
include_directories("${CMAKE_SOURCE_DIR}/src/core/misc")
include_directories("${CMAKE_SOURCE_DIR}/src/core/io")

add_executable(${PROJECT_NAME} test.cpp)

set(EXTRA_LIBS ${EXTRA_LIBS} cubic_spline)
set(EXTRA_LIBS ${EXTRA_LIBS} io)
set(EXTRA_LIBS ${EXTRA_LIBS} linalg)

target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS})
add_dependencies(${PROJECT_NAME} ${EXTRA_LIBS})

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
