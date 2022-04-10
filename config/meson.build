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

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel'
    lapack_vendor = 'mkl'
  endif
endif

if get_option('ilp64')
  add_project_arguments('-DWITH_ILP64', language:'fortran')
endif

if lapack_vendor == 'mkl'
  mkl_dep = []

  if get_option('ilp64')
    lpbit = 'i'
  else
    lpbit = ''
  endif

  if fc_id == 'intel'
    mkl_dep += fc.find_library('mkl_intel_'+lpbit+'lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += fc.find_library('mkl_gf_'+lpbit+'lp64')
    if get_option('openmp')
      mkl_dep += fc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += fc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += fc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'openblas'
  if get_option('ilp64')
    lpbit = '64'
  else
    lpbit = ''
  endif

  openblas_dep = dependency('openblas'+lpbit, required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library('openblas'+lpbit)
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack'+lpbit)
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  custom_libs = get_option('custom_libraries')
  if custom_libs[0].substring(0,2) == '-L'
    foreach lib: custom_libs
      if lib not in [custom_libs[0]]
        message(custom_libs[0].substring(2))
        lib_deps += fc.find_library(lib, dirs: custom_libs[0].substring(2))
      endif
    endforeach
  else
    foreach lib: custom_libs
      lib_deps += fc.find_library(lib)
    endforeach
  endif

else
  if get_option('ilp64')
    lpbit = '64'
  else
    lpbit = ''
  endif

  lapack_dep = dependency('lapack'+lpbit, required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack'+lpbit)
  endif
  lib_deps += lapack_dep
  blas_dep = dependency('blas'+lpbit, required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas'+lpbit)
  endif
  lib_deps += blas_dep
endif

# Create the tool chain library as subproject
mctc_dep = dependency('mctc-lib', required: false)
if not mctc_dep.found()
  mctc_prj = subproject(
    'mctc-lib',
    version: '>=0.2',
    default_options: [
      'default_library=static',
    ],
  )
  mctc_dep = mctc_prj.get_variable('mctc_dep')

  if install
    install_data(
      mctc_prj.get_variable('mctc_lic'),
      install_dir: get_option('datadir')/'licenses'/meson.project_name()/'mctc-lib'
    )
  endif
endif
lib_deps += mctc_dep
