# This file is part of multicharge. SPDX-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

cmake_minimum_required(VERSION 3.24)
get_directory_property(is-subproject PARENT_DIRECTORY)

project(
  "multicharge"
  LANGUAGES "Fortran"
  VERSION "0.5.0"
  DESCRIPTION "Molecular structure store for testing")

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Unprefixed options are kept for compatibility and only provide the defaults
option(
  WITH_OpenMP
  "Enable support for shared memory parallelisation with OpenMP (deprecated, use MULTICHARGE_WITH_OpenMP)"
  ON)
option(
  WITH_ILP64
  "Enable support for ILP64 BLAS/LAPACK calls (deprecated, use MULTICHARGE_WITH_ILP64)"
  OFF)
option(
  WITH_TESTS
  "Enable compilation of unit tests (deprecated, use MULTICHARGE_WITH_TESTS)"
  ON)
option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" OFF)
option(MULTICHARGE_WITH_OpenMP
       "Enable support for shared memory parallelisation with OpenMP"
       ${WITH_OpenMP})
option(MULTICHARGE_WITH_ILP64 "Enable support for ILP64 BLAS/LAPACK calls"
       ${WITH_ILP64})
option(MULTICHARGE_WITH_TESTS "Enable compilation of unit tests" ${WITH_TESTS})
option(MULTICHARGE_FETCH_DEPS
       "Download dependencies not found with find_package" ON)
set(MULTICHARGE_USE_MCTCLIB TRUE)

set("${PROJECT_NAME}-module-dir"
    "${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
    CACHE STRING "Subdirectory to install generated module files to")
set(module-dir "${${PROJECT_NAME}-module-dir}")
if(is-subproject)
  set(module-dir
      "${module-dir}"
      PARENT_SCOPE)
endif()

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Build type to be used." FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
endif()

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND MULTICHARGE_WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()

if(MULTICHARGE_WITH_ILP64)
  message(STATUS "Using LAPACK/BLAS ILP64 interface")
endif()

if(NOT TARGET "LAPACK::LAPACK")
  if(MULTICHARGE_WITH_ILP64)
    set(BLA_SIZEOF_INTEGER 8)
  endif()
  find_package("LAPACK" REQUIRED)
endif()

# Dependencies: use an installed package if available, otherwise download it
include(FetchContent)

set(deps_ "mctc-lib")

if(MULTICHARGE_WITH_TESTS)
  list(APPEND deps_ "mstore")
endif()

function(declare_dep name url tag)
  if(EXISTS "${PROJECT_SOURCE_DIR}/subprojects/${name}/CMakeLists.txt")
    FetchContent_Declare(
      "${name}" SOURCE_DIR "${PROJECT_SOURCE_DIR}/subprojects/${name}"
                FIND_PACKAGE_ARGS CONFIG)
  else()
    FetchContent_Declare(
      "${name}"
      GIT_REPOSITORY "${url}"
      GIT_TAG "${tag}"
      GIT_SHALLOW TRUE
      FIND_PACKAGE_ARGS CONFIG)
  endif()
endfunction()

declare_dep("mctc-lib" "https://github.com/grimme-lab/mctc-lib" "v0.6.1")

if(MULTICHARGE_WITH_TESTS)
  declare_dep("mstore" "https://github.com/grimme-lab/mstore" "v0.4.0")
endif()

FetchContent_MakeAvailable(${deps_})

unset(deps_)

# mstore does not define a namespaced target when built in-tree
if(TARGET "mstore" AND NOT TARGET "mstore::mstore")
  add_library("mstore::mstore" ALIAS "mstore")
endif()

# Collect source of the project
add_subdirectory("src")

# Add example application
add_subdirectory("app")

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
add_library("${PROJECT_NAME}::${PROJECT_NAME}" ALIAS "${PROJECT_NAME}")
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")

install(TARGETS "${PROJECT_NAME}" "${PROJECT_NAME}-lib"
        EXPORT "${PROJECT_NAME}-targets")

install(
  EXPORT "${PROJECT_NAME}-targets"
  NAMESPACE "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}")

# Package license files
install(FILES "LICENSE"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}")

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config/template-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config/template.pc"
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

# add the testsuite
if(MULTICHARGE_WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
endif()
