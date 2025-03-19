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

if ((MULTICHARGE_BLAS MATCHES ^Intel) OR (DEFINED ENV{MKLROOT}))
  enable_language("C")
endif()

if(WITH_ILP64)
  set(BLA_SIZEOF_INTEGER 8)
  set(_nvpl_int "_ilp64")
else()
  set(_nvpl_int "_lp64")
  set(BLA_SIZEOF_INTEGER 4)
endif()

if(NOT MULTICHARGE_BLAS_FOUND)
  if(MULTICHARGE_BLAS STREQUAL "NVPL")
    find_package("nvpl_blas" REQUIRED)
    set(MULTICHARGE_BLAS_FOUND TRUE)

    if((BLA_THREAD STREQUAL "OMP") OR (BLA_THREAD STREQUAL "ANY"))
      set(_nvpl_thread "_omp")
    else()
      set(_nvpl_thread "_seq")
    endif()

    add_library("multicharge::BLAS" INTERFACE IMPORTED GLOBAL)
    target_link_libraries("multicharge::BLAS" INTERFACE "nvpl::blas${_nvpl_int}${_nvpl_thread}")
  else()
    find_package("BLAS" REQUIRED)
    set(MULTICHARGE_BLAS_FOUND ${BLAS_FOUND})
  
    if(NOT TARGET "multicharge::BLAS")
      add_library("multicharge::BLAS" INTERFACE IMPORTED GLOBAL)
      target_link_libraries("multicharge::BLAS" INTERFACE "BLAS::BLAS")
    endif()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(multicharge-blas DEFAULT_MSG MULTICHARGE_BLAS_FOUND)
