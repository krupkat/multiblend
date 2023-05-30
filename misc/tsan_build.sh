cmake -B build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_TOOLCHAIN_FILE="vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_TARGET_TRIPLET="x64-linux-tsan" \
  -DMULTIBLEND_WITH_TSAN=ON

cmake --build build -j $(nproc)

ctest --test-dir build -V
