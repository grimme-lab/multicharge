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

# Include the mctc-lib project
if(NOT TARGET mctc-lib)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mctc-lib/CMakeLists.txt)
    add_subdirectory("mctc-lib")
  else()
    set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
    message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mctc-lib"
      GIT_REPOSITORY "${mctc-lib-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mctc-lib")
  endif()
endif()

list(
  APPEND lib-deps
  "mctc-lib"
)

# Include the mstore project
if(NOT TARGET mstore)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mstore/CMakeLists.txt)
    add_subdirectory("mstore")
  else()
    set("mstore-url" "https://github.com/grimme-lab/mstore")
    message(STATUS "Retrieving mstore from ${mstore-url}")
    include(FetchContent)
    FetchContent_Declare(
      "mstore"
      GIT_REPOSITORY "${mstore-url}"
      GIT_TAG "HEAD"
    )
    FetchContent_MakeAvailable("mstore")
  endif()
endif()

set(lib-deps "${lib-deps}" PARENT_SCOPE)
