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

validation_tester = find_program(files('tester.py'))

run_element = [
  '01_h3+',
  '02_heada',
  '03_etlicl-',
]

foreach t : run_element
  benchmark(
    '@0@'.format(t),
    validation_tester,
    args: [
      multicharge_exe,
      files(t/'multicharge_ref.json'),
      files(t/'coord'),
      '--grad',
      t.contains('+') ? ['--charge', '+1'] : [],
      t.contains('-') ? ['--charge', '-1'] : [],
    ],
    suite: 'multicharge',
  )
endforeach
