@PACKAGE_INIT@

set(MULTICHARGE_WITH_OpenMP @MULTICHARGE_WITH_OpenMP@)
set(MULTICHARGE_WITH_ILP64 @MULTICHARGE_WITH_ILP64@)
set(MULTICHARGE_USE_MCTCLIB @MULTICHARGE_USE_MCTCLIB@)

enable_language("Fortran")

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
  include(CMakeFindDependencyMacro)

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND MULTICHARGE_WITH_OpenMP)
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "LAPACK::LAPACK")
    find_dependency("LAPACK")
  endif()

  if(NOT TARGET "mctc-lib::mctc-lib" AND MULTICHARGE_USE_MCTCLIB)
    find_dependency("mctc-lib")
  endif()

endif()
