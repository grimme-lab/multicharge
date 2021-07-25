# This file is part of multicharge.
# SPDX-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Handling of subproject dependencies
macro(
  "multicharge_find_package"
  package
  methods
  url
)
  string(TOLOWER "${package}" _pkg_lc)
  string(TOUPPER "${package}" _pkg_uc)

  foreach(method in ITEMS ${methods})

    if("${method}" STREQUAL "cmake")
      message(STATUS "${package}: Find installed package")
      find_package("${package}" CONFIG)
      if("${package}_FOUND")
        message(STATUS "${package}: Found installed package")
        break()
      endif()
    endif()

    if("${method}" STREQUAL "pkgconf")
      find_package(PkgConfig QUIET)
      pkg_check_modules("${_pkg_uc}" QUIET "${package}")
      if("${_pkg_uc}_FOUND")
        message(STATUS "Found ${package} via pkg-config")

        add_library("${package}::${package}" INTERFACE IMPORTED)
        target_link_libraries(
          "${package}::${package}"
          INTERFACE
          "${${_pkg_uc}_LINK_LIBRARIES}"
        )
        target_include_directories(
          "${package}::${package}"
          INTERFACE
          "${${_pkg_uc}_INCLUDE_DIRS}"
        )
        break()
      endif()
    endif()

    if("${method}" STREQUAL "subproject")
      set("${_pkg_uc}_SOURCE_DIR" "${PROJECT_SOURCE_DIR}/subprojects/${package}")
      set("${_pkg_uc}_BINARY_DIR" "${PROJECT_BINARY_DIR}/subprojects/${package}")
      if(EXISTS "${${_pkg_uc}_SOURCE_DIR}/CMakeLists.txt")
        message(STATUS "Include ${package} from subprojects")
        add_subdirectory(
          "${${_pkg_uc}_SOURCE_DIR}"
          "${${_pkg_uc}_BINARY_DIR}"
        )

        add_library("${package}::${package}" INTERFACE IMPORTED)
        target_link_libraries("${package}::${package}" INTERFACE "${package}")

        # We need the module directory in the subproject before we finish the configure stage
        if(NOT EXISTS "${${_pkg_uc}_BINARY_DIR}/include")
          make_directory("${${_pkg_uc}_BINARY_DIR}/include")
        endif()

        break()
      endif()
    endif()

    if("${method}" STREQUAL "fetch")
      message(STATUS "Retrieving ${package} from ${url}")
      include(FetchContent)
      FetchContent_Declare(
        "${_pkg_lc}"
        GIT_REPOSITORY "${url}"
        GIT_TAG "HEAD"
      )
      FetchContent_MakeAvailable("${_pkg_lc}")

      add_library("${package}::${package}" INTERFACE IMPORTED)
      target_link_libraries("${package}::${package}" INTERFACE "${package}")

      # We need the module directory in the subproject before we finish the configure stage
      FetchContent_GetProperties("${_pkg_lc}" BINARY_DIR "${_pkg_uc}_BINARY_DIR")
      if(NOT EXISTS "${${_pkg_uc}_BINARY_DIR}/include")
        make_directory("${${_pkg_uc}_BINARY_DIR}/include")
      endif()

      break()
    endif()

  endforeach()

  unset(_pkg_lc)
  unset(_pkg_uc)

  if(NOT TARGET "${package}::${package}")
    message(FATAL_ERROR "Could not find dependency ${package}")
  endif()

endmacro()
