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

os = host_machine.system()
fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

blas_dep = dependency(
  'blas-interface',
  fallback: ['blas-interface', 'blas_interface_dep'],
  default_options: [
    'default_library=static',
  ],
)
lib_deps += blas_dep

lapack_dep = dependency(
  'lapack-interface',
  fallback: ['lapack-interface', 'lapack_interface_dep'],
  default_options: [
    'default_library=static',
  ],
)
lib_deps += lapack_dep

# Create the tool chain library
mctc_dep = dependency(
  'mctc-lib',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: [
    'default_library=static',
  ],
)
lib_deps += mctc_dep
