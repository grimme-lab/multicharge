@PACKAGE_INIT@

set("@PROJECT_NAME@_WITH_OpenMP" @WITH_OpenMP@)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")

  include(CMakeFindDependencyMacro)

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND "@PROJECT_NAME@_WITH_OpenMP")
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "LAPACK::LAPACK")
    find_dependency("LAPACK")
  endif()

  if(NOT TARGET "mctc-lib::mctc-lib")
    find_dependency("mctc-lib")
  endif()
endif()
