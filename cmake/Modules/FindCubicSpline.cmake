find_path(CubicSpline_INCLUDE_DIR CubicSpline.h
          HINTS contrib
          PATH_SUFFIXES CubicSpline)

# find_library(LIBXML2_LIBRARY NAMES xml2 libxml2
#              HINTS ${PC_LIBXML_LIBDIR} ${PC_LIBXML_LIBRARY_DIRS} )

# include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
# find_package_handle_standard_args(CubicSpline  DEFAULT_MSG
#                                   CubicSpline_LIBRARY CubicSpline_INCLUDE_DIR)

# mark_as_advanced(CubicSpline_INCLUDE_DIR CubicSpline_LIBRARY)

# set(CubicSpline_LIBRARIES ${CubicSpline_LIBRARY})
set(CubicSpline_INCLUDE_DIRS ${CubicSpline_INCLUDE_DIR})