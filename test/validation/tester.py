#!/usr/bin/env python

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
"""
Minimal Python wrapper for testing the multicharge command line interface.

The wrapper will assume a specific order in the arguments rather than
providing a generic command line interface by itself since it is
supposed to be used by meson for testing purposes only.
"""

try:
    import json
    import os
    import subprocess
    import sys

    import numpy as np
    import pytest
except ImportError:
    exit(77)

if len(sys.argv) < 4:
    raise RuntimeError("Requires at least four arguments")

THR = 5.0e-7
prog = sys.argv[1]
reference_file = sys.argv[2]
args = sys.argv[3:]

method_name = os.path.basename(reference_file)
test_name = os.path.basename(os.path.dirname(reference_file))

stat = subprocess.call(
    [prog] + args + ["--json"],
    shell=False,
    stdin=None,
    stderr=subprocess.STDOUT,
)
if stat != 0:
    raise RuntimeError("Calculation failed")

with open(reference_file, encoding="utf8") as f:
    ref = json.load(f)
    del ref["version"]

with open("multicharge.json", encoding="utf8") as f:
    res = json.load(f)

for key in ref:
    if key not in res:
        raise RuntimeError("Missing '" + key + "' entry in results")
    _res = np.array(res[key])
    _ref = np.array(ref[key])
    assert pytest.approx(_res, abs=THR) == _ref, \
        f"mismatch for {key}:\n" \
        f"+actual\n{_res}\n" \
        f"-reference\n{_ref}\n" \
        f"@difference\n{_res - _ref}"
