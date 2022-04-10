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

subdir('data')
subdir('param')

if get_option('ilp64')
  lpbit = '64'
else
  lpbit = ''
endif

srcs += files(
  'blas'+lpbit+'.f90',
  'cutoff.f90',
  'data.f90',
  'ewald.f90',
  'lapack'+lpbit+'.f90',
  'model.F90',
  'ncoord.f90',
  'output.f90',
  'param.f90',
  'version.f90',
  'wignerseitz.f90',
)
