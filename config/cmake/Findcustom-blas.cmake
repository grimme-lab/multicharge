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
endif()

if(NOT BLAS_FOUND)
  find_package("BLAS")
  if(NOT TARGET "BLAS::BLAS")
    add_library("BLAS::BLAS" INTERFACE IMPORTED)
    target_link_libraries("BLAS::BLAS" INTERFACE "${BLAS_LIBRARIES}")
  endif()
endif()
