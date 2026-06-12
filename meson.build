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

project(
  'multicharge',
  'fortran',
  version: '0.5.0',
  license: 'Apache-2.0',
  meson_version: '>=0.60',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: false, native: false)

# Check for specific unsupported meson versions
if meson.version().version_compare('==1.8.0')
  error('Meson version 1.8.0 has a known issue — please use any other version ≥ 0.55.0')
endif

# General configuration information
inc_dirs = []
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
multicharge_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
multicharge_inc = multicharge_lib.private_dir_include()
multicharge_dep = declare_dependency(
  link_with: multicharge_lib,
  include_directories: multicharge_inc,
  dependencies: lib_deps,
)

# Add applications
subdir('app')

# Package the license files
multicharge_lic = files(
  'LICENSE',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    multicharge_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    multicharge_lib,
    description: 'Electronegativity equilibration model for charges and multipoles',
    subdirs: ['', module_id],
  )

  cmake_data = configuration_data({
    'PROJECT_NAME': meson.project_name(),
    'PROJECT_VERSION': meson.project_version(),
    'PACKAGE_INIT': '',
    'CMAKE_INSTALL_PREFIX': get_option('prefix'),
    'CMAKE_INSTALL_INCLUDEDIR': get_option('includedir'),
    'CMAKE_INSTALL_LIBDIR': get_option('libdir'),
    'MULTICHARGE_WITH_OpenMP': get_option('openmp') ? 'ON' : 'OFF',
    'MULTICHARGE_WITH_ILP64': get_option('ilp64') ? 'ON' : 'OFF',
    'MULTICHARGE_USE_MCTCLIB': mctc_dep.type_name() != 'internal' ? 'ON' : 'OFF',
    'module-dir': module_id,
  })
  cmake_files = [
    configure_file(
      input: files('config'/'template-config.cmake'),
      output: meson.project_name() + '-config.cmake',
      configuration: cmake_data,
    ),
    configure_file(
      input: files('config'/'template-targets.cmake'),
      output: meson.project_name() + '-targets.cmake',
      configuration: cmake_data,
    ),
    files('config'/'cmake'/'Findcustom-blas.cmake'),
    files('config'/'cmake'/'Findcustom-lapack.cmake'),
    files('config'/'cmake'/'Findmctc-lib.cmake'),
    files('config'/'cmake'/'Findmstore.cmake'),
    files('config'/'cmake'/'multicharge-utils.cmake'),
  ]
  install_data(
    cmake_files,
    install_dir: get_option('libdir')/'cmake'/meson.project_name(),
    install_tag: 'devel',
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man/multicharge.1.adoc'),
        output: '@BASENAME@',
      )
    )
  endif
endif

# add the testsuite
subdir('test')
