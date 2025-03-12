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

if ((BLA_VENDOR MATCHES ^Intel) OR (DEFINED ENV{MKLROOT}))
  enable_language("C")
endif()

if(WITH_ILP64)
  set(BLA_SIZEOF_INTEGER 8)
  set(_nvpl_int "_ilp64")
else()
  set(_nvpl_int "_lp64")
endif()

if(NOT LAPACK_FOUND)
  if(BLA_VENDOR STREQUAL "NVPL")
    find_package("nvpl_lapack" REQUIRED)

    if((BLA_THREAD STREQUAL "OMP") OR (BLA_THREAD STREQUAL "ANY"))
      set(_nvpl_thread "_omp")
    else()
      set(_nvpl_thread "_seq")
    endif()

    add_library("LAPACK::LAPACK" INTERFACE IMPORTED)
    target_link_libraries("LAPACK::LAPACK" INTERFACE "nvpl::lapack${_nvpl_int}${_nvpl_thread}")
  else()
    find_package("LAPACK" REQUIRED)
  endif()

  if(NOT TARGET "BLAS::BLAS")
    find_package("custom-blas" REQUIRED)
  endif()

  if(NOT TARGET "LAPACK::LAPACK")
    add_library("LAPACK::LAPACK" INTERFACE IMPORTED)
    target_link_libraries("LAPACK::LAPACK" INTERFACE "${LAPACK_LIBRARIES}" "BLAS::BLAS")
  endif()
endif()
