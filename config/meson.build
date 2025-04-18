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
cc = fc
if has_cc
  cc = meson.get_compiler('c')
endif
fc_id = fc.get_id()

# Treat Intel LLVM compiler as standard Intel compiler for compatibility
if fc_id == 'intel-llvm'
  fc_id = 'intel'
endif
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
  ilp64 = true
  add_project_arguments('-DIK=i8', language:'fortran')
else
  add_project_arguments('-DIK=i4', language:'fortran')
  ilp64 = false
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel'
    mkl_dep += cc.find_library(ilp64 ? 'mkl_intel_ilp64' : 'mkl_intel_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += cc.find_library(ilp64 ? 'mkl_gf_ilp64' : 'mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += cc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += cc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = cc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'openblas'
  openblas_dep = dependency(ilp64 ? 'openblas64' : 'openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = cc.find_library(ilp64 ? 'openblas64' : 'openblas') 
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency(ilp64 ? 'lapack64' : 'lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = cc.find_library(ilp64 ? 'lapack64' : 'lapack')
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  custom_deps = []
  libs = get_option('custom_libraries')
  if libs.length() == 0
    error('Custom LAPACK vendor selected, but no custom libraries were provided. Please specify the custom libraries using the "-Dcustom_libraries" option.')
  endif
  if libs[0].startswith('-L')
    foreach lib: libs
      if lib != libs[0]
        custom_deps += cc.find_library(lib, dirs: libs[0].substring(2))
      endif
    endforeach
  else
    foreach lib: libs
      custom_deps += cc.find_library(lib)
    endforeach
  endif
  if (not fc.links('external dsytrs; call dsytrs(); end', dependencies: (get_option('openmp') ? [custom_deps, omp_dep] : [custom_deps])))
    error('Custom LAPACK libraries do not link')
  elif (not fc.links('external dsytrs; call dgemm(); end', dependencies: (get_option('openmp') ? [custom_deps, omp_dep] : [custom_deps])))
    error('Custom BLAS libraries do not link')
  endif
  lib_deps += custom_deps


else
  lapack_dep = dependency(ilp64 ? 'lapack64' : 'lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = cc.find_library(ilp64 ? 'lapack64' : 'lapack') 
  endif
  lib_deps += lapack_dep
  blas_dep = dependency(ilp64 ? 'blas64' : 'blas', required: false) 
  if not blas_dep.found()
    blas_dep = cc.find_library(ilp64 ? 'blas64' : 'blas') 
  endif
  lib_deps += blas_dep
endif

# Create the tool chain library as subproject
mctc_dep = dependency('mctc-lib', required: false)
if not mctc_dep.found()
  mctc_prj = subproject(
    'mctc-lib',
    version: '>=0.4.1',
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
