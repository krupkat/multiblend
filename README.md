[![tests](https://github.com/krupkat/multiblend/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/krupkat/multiblend/actions/workflows/build.yml)
[![clang-format](https://github.com/krupkat/multiblend/actions/workflows/clang-format-check.yml/badge.svg?branch=main)](https://github.com/krupkat/multiblend/actions/workflows/clang-format-check.yml)
[![clang-tidy](https://github.com/krupkat/multiblend/actions/workflows/clang-tidy-check.yml/badge.svg?branch=main)](https://github.com/krupkat/multiblend/actions/workflows/clang-tidy-check.yml)

# Multiblend
Fork of multiblend: http://horman.net/multiblend/

Multiblend 2.0 (c) 2021 David Horman

With modifications (c) 2023 Tomas Krupka

## Changes

Changes from version 2.0rc5:

- Added CMake + isolated the algorithm in a cmake target
  - Translation units: https://github.com/krupkat/multiblend/pull/3
  - Separate algorithm: https://github.com/krupkat/multiblend/pull/4
  - Cmake target: https://github.com/krupkat/multiblend/pull/17
- Added CI with tests, clang-format and clang-tidy
  - Run test in CI: https://github.com/krupkat/multiblend/pull/1
  - clang-format: https://github.com/krupkat/multiblend/pull/2
  - clang-tidy: https://github.com/krupkat/multiblend/pull/7
  - Tests: https://github.com/krupkat/multiblend/tree/main/tests
- Fixed memory leaks
  - General + jpg: https://github.com/krupkat/multiblend/pull/10
  - Png: https://github.com/krupkat/multiblend/pull/14
  - Tiff: https://github.com/krupkat/multiblend/pull/15
  - new + malloc: https://github.com/krupkat/multiblend/pull/16
- Replaced logging with [spdlog](https://github.com/gabime/spdlog/)
  - Spdlog library + exceptions: https://github.com/krupkat/multiblend/pull/18
- Replaced multi threading with [threadpool](https://github.com/bshoshany/thread-pool)
  - Threadpool library: https://github.com/krupkat/multiblend/pull/19
  - Remove global threadpool: https://github.com/krupkat/multiblend/pull/22
- Further cleanup
  - Namespaces: https://github.com/krupkat/multiblend/pull/5
  - Member naming: https://github.com/krupkat/multiblend/pull/6
  - Remove macros: https://github.com/krupkat/multiblend/pull/8
  - Disable MapAlloc: https://github.com/krupkat/multiblend/pull/20
  - C++ memory alignment: https://github.com/krupkat/multiblend/pull/21
- Non x86 support
  - SIMDe integration: https://github.com/krupkat/multiblend/pull/25

## Build

Requires CMake 3.21 and a C++20 compiler. Use the included `sh` scripts on Linux and the `ps1` scripts on Windows. First download [vcpkg](https://github.com/microsoft/vcpkg):

```
./misc/get_vcpkg.sh
```

Regular build:

```
./misc/build.sh
```

Linux / Mac only:

Build with [address sanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer):

```
./misc/asan_build.sh
```

Build with [thread sanitizer](https://github.com/google/sanitizers/wiki/ThreadSanitizerCppManual):

```
./misc/tsan_build.sh
```

## Test requirements

Python + pip

## Library integration

Check out [Xpano](https://github.com/krupkat/xpano) for an example of integrating the Multiblend library in your C++ project.
