cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_TARGET_TRIPLET="x64-linux-asan" \
  -DMULTIBLEND_WITH_ASAN=ON

cmake --build build -j $(nproc)

ctest --test-dir build -V
