project(io)

include_directories("${CMAKE_SOURCE_DIR}/src/core/misc")
add_library(${PROJECT_NAME} STATIC io.cpp)
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
