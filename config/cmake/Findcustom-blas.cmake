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

if(NOT BLAS_FOUND)

  if (WITH_ILP64)
    if(CMAKE_VERSION VERSION_LESS "3.13")
      if(NOT LAPACK_LIBRARIES)
        message(FATAL_ERROR "Use of WITH_ILP64 with CMAKE_VERSION < 3.13 requires user to set LAPACK_LIBRARIES")
      endif()
    elseif(CMAKE_VERSION VERSION_LESS "3.22")
      if (NOT LAPACK_LIBRARIES)
        if(NOT ((BLA_VENDOR MATCHES "64ilp") OR (BLA_VENDOR MATCHES "_dyn")))
          message(FATAL_ERROR
            "Use of WITH_ILP64 with CMAKE_VERSION < 3.22 requires user to set BLA_VENDOR to a valid ilp64 provider or to set LAPACK_LIBRARIES")
        endif()
      endif()
    endif()
  endif()
        
  find_package("BLAS")

  if(NOT TARGET "BLAS::BLAS")
    add_library("BLAS::BLAS" INTERFACE IMPORTED)
    target_link_libraries("BLAS::BLAS" INTERFACE "${BLAS_LIBRARIES}")
  endif()

endif()
