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

set(srcs)
add_subdirectory("multicharge")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")
list(
  APPEND srcs
  "${dir}/multicharge.f90"
)

# multicharge library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
add_library("${PROJECT_NAME}::${PROJECT_NAME}-lib" ALIAS "${PROJECT_NAME}-lib")

set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)

target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "mctc-lib::mctc-lib"
  "LAPACK::LAPACK"
  "$<$<BOOL:${MULTICHARGE_WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
)

target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "IK=$<IF:$<BOOL:${MULTICHARGE_WITH_ILP64}>,i8,i4>"
)

target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
