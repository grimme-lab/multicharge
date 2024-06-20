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

add_subdirectory("data")
add_subdirectory("param")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/blas.f90"
  "${dir}/cutoff.f90"
  "${dir}/data.f90"
  "${dir}/ewald.f90"
  "${dir}/lapack.f90"
  "${dir}/model.f90"
  "${dir}/ncoord.f90"
  "${dir}/output.f90"
  "${dir}/param.f90"
  "${dir}/version.f90"
  "${dir}/wignerseitz.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
