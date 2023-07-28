cmake -B build -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_TOOLCHAIN_FILE="vcpkg/scripts/buildsystems/vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET="x64-windows-static" `
  -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded"

cmake --build build

ctest --test-dir build -V
