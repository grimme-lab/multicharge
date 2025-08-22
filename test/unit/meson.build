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

# Find mstore dependency for testing
mstore_dep = dependency(
  'mstore',
  version: '>=0.3.0',
  fallback: ['mstore', 'mstore_dep'],
  required: not meson.is_subproject(),
  default_options: [
    'default_library=static',
  ],
)
# If we do not find mstore and are a subproject, we just skip testing
if not mstore_dep.found()
   subdir_done()
endif

tests = [
  'model',
  'pbc',
  'wignerseitz',
]

test_srcs = files(
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [multicharge_dep, mstore_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(t, tester, args: t)
endforeach
