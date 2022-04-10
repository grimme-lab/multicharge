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

function(check_lapack)
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LINK_OPTIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
check_fortran_source_compiles("
      program check_lapack
      implicit none (type, external)
      external :: disnan
      end program" LAPACK_links)
set(LAPACK_links ${LAPACK_links} PARENT_SCOPE)
endfunction(check_lapack)

function(check_blas)
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LINK_OPTIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
check_fortran_source_compiles("
      program check_blas
      implicit none (type, external)
      external :: dgemm
      end program" BLAS_links)
set(BLAS_links ${BLAS_links} PARENT_SCOPE)
endfunction(check_blas)


if(NOT LAPACK_FOUND)

  # Safeguard in case a vendor was specified
  if(BLA_VENDOR MATCHES "_64ilp")
    if(NOT WITH_ILP64)
      message(STATUS "BLA_VENDOR ${BLA_VENDOR} uses 64-bit integers. Setting WITH_ILP64")
      set(WITH_ILP64 ON PARENT_SCOPE)
    endif()
  endif()

  # ILP64 check
  if (WITH_ILP64 AND (NOT LAPACK_LIBRARIES))
    if(NOT ((BLA_VENDOR MATCHES "_64ilp") OR (BLA_VENDOR MATCHES "_dyn")))
      if(CMAKE_VERSION VERSION_LESS "3.22")
        message(FATAL_ERROR
          "Use of WITH_ILP64 with CMAKE_VERSION < 3.22 requires user to set BLA_VENDOR to a valid ilp64 provider or to set LAPACK_LIBRARIES")
      elseif(NOT BLA_SIZEOF_INTEGER EQUAL 8)
        set(BLA_SIZEOF_INTEGER 8 FORCE)
      endif()
    endif()
  endif()

  if(LAPACK_LIBRARIES)

    # Check if LAPACK links with the passed libraries
    if(NOT TARGET LAPACK::LAPACK)
      check_lapack()
      if(LAPACK_links)
        set (LAPACK_links ${LAPACK_links})
        add_library(LAPACK::LAPACK INTERFACE IMPORTED)
        set_target_properties(LAPACK::LAPACK PROPERTIES
          INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
      else()
        message(FATAL_ERROR "Custom LAPACK library does not link")
      endif()
    endif()

    # Check if BLAS links with the same set of libraries
    if(NOT TARGET BLAS::BLAS)
      check_blas()
      if(BLAS_links)
        set (BLAS_links ${BLAS_links})
        add_library(BLAS::BLAS INTERFACE IMPORTED)
      else()
        message(FATAL_ERROR "Custom LAPACK library must contain BLAS")
      endif()    
    endif()

    set(LAPACK_FOUND TRUE)

  else()

    find_package("LAPACK")

    if(NOT TARGET "BLAS::BLAS")
      find_package("custom-blas")
    endif()

    if(NOT TARGET "LAPACK::LAPACK")
      add_library("LAPACK::LAPACK" INTERFACE IMPORTED)
      target_link_libraries("LAPACK::LAPACK" INTERFACE "${LAPACK_LIBRARIES}" "BLAS::BLAS")
    endif()

  endif()
endif()
