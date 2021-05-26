# Atomic Partial Charge Model

[![Apache-2.0](https://img.shields.io/github/license/grimme-lab/multicharge)](LICENSE)
[![Release](https://img.shields.io/github/v/release/grimme-lab/multicharge)](https://github.com/grimme-lab/multicharge/releases/latest)
[![CI](https://github.com/grimme-lab/multicharge/workflows/CI/badge.svg)](https://github.com/grimme-lab/multicharge/actions)
[![codecov](https://codecov.io/gh/grimme-lab/multicharge/branch/main/graph/badge.svg)](https://codecov.io/gh/grimme-lab/multicharge)

Electronegativity equilibration model for atomic partial charges.


## Installation

To build this project from the source code in this repository you need to have
a Fortran compiler supporting Fortran 2008 and one of the supported build systems:
- [meson](https://mesonbuild.com) version 0.53 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
- [cmake](https://cmake.org) version 3.14 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- [fpm](https://github.com/fortran-lang/fpm) version 0.2.0 or newer

Currently this project supports GCC and Intel compilers.


### Building with meson

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile the project run

```
meson compile -C _build
```

You can run the projects testsuite with

```
meson test -C _build --print-errorlogs
```

To include ``multicharge`` in your project add the following wrap file to your subprojects directory:

```ini
[wrap-git]
directory = multicharge
url = https://github.com/grimme-lab/multicharge
revision = head
```

You can retrieve the dependency from the wrap fallback with

```meson
multicharge_dep = dependency('multicharge', ['multicharge', 'multicharge_dep'])
```

and add it as dependency to your targets.


### Building with CMake

Alternatively, this project can be build with CMake (in this case ninja 1.10 or newer is required):

```
cmake -B _build -G Ninja
```

To compile the project with CMake run

```
cmake --build _build
```

You can run the project testsuite with

```
pushd _build && ctest && popd
```

To include ``multicharge`` in your CMake project retrieve it using the ``FetchContent`` module:

```cmake
if(NOT TARGET multicharge)
  set("multicharge-url" "https://github.com/grimme-lab/multicharge")
  message(STATUS "Retrieving multicharge from ${multicharge-url}")
  include(FetchContent)
  FetchContent_Declare(
    "multicharge"
    GIT_REPOSITORY "${multicharge-url}"
    GIT_TAG "HEAD"
  )
  FetchContent_MakeAvailable("multicharge")
endif()
```

And link against the ``"multicharge"`` interface library.

```cmake
target_link_libraries("${PROJECT_NAME}-lib" PUBLIC "multicharge")
```


### Building with fpm

Invoke fpm in the project root with

```
fpm build
```

To run the testsuite use

```
fpm test
```

You can access the ``multicharge`` program using the run subcommand

```
fpm run -- --help
```

To use ``multicharge`` for testing include it as dependency in your package manifest

```toml
[dependencies]
multicharge.git = "https://github.com/grimme-lab/multicharge"
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
