cmake -B build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_TOOLCHAIN_FILE="vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_TARGET_TRIPLET="arm64-osx" \
  -DMULTIBLEND_WITH_NEON=ON

cmake --build build -j 6

ctest --test-dir build -V
