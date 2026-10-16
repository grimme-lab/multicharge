# Imported targets for a meson installation, the CMake build exports its own targets.
# PACKAGE_PREFIX_DIR is set relative to this file by the package config file.

set(CMAKE_IMPORT_FILE_VERSION 1)

find_library(
  MULTICHARGE_LIBRARY
  "@PROJECT_NAME@"
  PATHS "${PACKAGE_PREFIX_DIR}/@CMAKE_INSTALL_LIBDIR@"
  NO_DEFAULT_PATH
)
if(NOT MULTICHARGE_LIBRARY)
  message(FATAL_ERROR "@PROJECT_NAME@ library not found in ${PACKAGE_PREFIX_DIR}/@CMAKE_INSTALL_LIBDIR@")
endif()

set(_multicharge_deps "LAPACK::LAPACK")
if(MULTICHARGE_USE_MCTCLIB)
  list(APPEND _multicharge_deps "mctc-lib::mctc-lib")
endif()
if(MULTICHARGE_WITH_OpenMP)
  list(APPEND _multicharge_deps "OpenMP::OpenMP_Fortran")
endif()

# Shared or static depends on what meson installed and what find_library prefers
add_library("@PROJECT_NAME@::@PROJECT_NAME@-lib" UNKNOWN IMPORTED)
set_target_properties(
  "@PROJECT_NAME@::@PROJECT_NAME@-lib"
  PROPERTIES
  IMPORTED_LOCATION "${MULTICHARGE_LIBRARY}"
  IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran"
  INTERFACE_INCLUDE_DIRECTORIES
  "${PACKAGE_PREFIX_DIR}/@CMAKE_INSTALL_INCLUDEDIR@;${PACKAGE_PREFIX_DIR}/@CMAKE_INSTALL_INCLUDEDIR@/@module-dir@"
  INTERFACE_LINK_LIBRARIES "${_multicharge_deps}"
)
unset(_multicharge_deps)

add_library("@PROJECT_NAME@::@PROJECT_NAME@" INTERFACE IMPORTED)
target_link_libraries("@PROJECT_NAME@::@PROJECT_NAME@" INTERFACE "@PROJECT_NAME@::@PROJECT_NAME@-lib")

set(CMAKE_IMPORT_FILE_VERSION)
