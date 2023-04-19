cmake_minimum_required(VERSION 3.21)

project(Multiblend)
include(CTest)

set(CMAKE_CXX_STANDARD 20)

add_executable(Multiblend
  "src/multiblend.cpp"
)

find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)
find_package(TIFF REQUIRED)

target_include_directories(Multiblend PRIVATE 
  "src"
)

target_link_libraries(Multiblend
  JPEG::JPEG
  PNG::PNG
  TIFF::TIFF
)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
