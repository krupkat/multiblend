export VCPKG_FORCE_SYSTEM_BINARIES=1 # needed for vcpkg on arm

cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_TARGET_TRIPLET="arm64-osx" \
  -DMULTIBLEND_ARM_OPTIMIZED=ON

cmake --build build -j 6

ctest --test-dir build -V
