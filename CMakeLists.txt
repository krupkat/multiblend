cmake_minimum_required(VERSION 3.21)

project(Multiblend)
include(CTest)

set(CMAKE_CXX_STANDARD 20)

add_executable(Multiblend
  "src/functions.cpp"
  "src/geotiff.cpp"
  "src/image.cpp"
  "src/main.cpp"
  "src/mapalloc.cpp"
  "src/multiblend.cpp"
  "src/pnger.cpp"
  "src/pyramid.cpp"
  "src/threadpool.cpp"
)

find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)
find_package(TIFF REQUIRED)

target_include_directories(Multiblend PRIVATE 
  "."
)

target_link_libraries(Multiblend
  JPEG::JPEG
  PNG::PNG
  TIFF::TIFF
)

if(NOT MSVC)
    target_compile_options(Multiblend PRIVATE -msse4.1)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
