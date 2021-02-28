# Atomic Partial Charge Model

[![Apache-2.0](https://img.shields.io/github/license/grimme-lab/multicharge)](LICENSE)
[![Release](https://img.shields.io/github/v/release/grimme-lab/multicharge)](https://github.com/grimme-lab/multicharge/releases/latest)
[![CI](https://github.com/grimme-lab/multicharge/workflows/CI/badge.svg)](https://github.com/grimme-lab/multicharge/actions)
[![codecov](https://codecov.io/gh/grimme-lab/multicharge/branch/main/graph/badge.svg)](https://codecov.io/gh/grimme-lab/multicharge)

Electronegativity equilibration model for atomic partial charges.


## Installation

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.53 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable, currently this project supports GCC and Intel compilers.
To compile the project run

```
meson compile -C _build
```

You can run the projects testsuite with

```
meson test -C _build --print-errorlogs
```


## Usage

For an overview over all command line arguments use the ``--help`` argument or checkout the [``multicharge(1)``](man/multicharge.1.adoc) manpage.


## License

Licensed under the Apache License, Version 2.0 (the “License”);
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an *“as is” basis*,
*without warranties or conditions of any kind*, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Apache-2.0 license, shall be licensed as above, without any additional
terms or conditions.
