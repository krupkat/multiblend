cmake_minimum_required(VERSION 3.21)
project(Multiblend)

set(MULTIBLEND_ROOT_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MULTIBLEND_ROOT_PROJECT ON)
endif()

OPTION(MULTIBLEND_WITH_ASAN "Build with address sanitizer" OFF)
OPTION(MULTIBLEND_BUILD_EXE "Build executable" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_LIB_WITH_JPEG "Build with JPEG support" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_LIB_WITH_PNG "Build with PNG support" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_LIB_WITH_TIFF "Build with TIFF support" ${MULTIBLEND_ROOT_PROJECT})

if(MULTIBLEND_BUILD_EXE)
  if(NOT MULTIBLEND_LIB_WITH_JPEG OR NOT MULTIBLEND_LIB_WITH_PNG OR NOT MULTIBLEND_LIB_WITH_TIFF)
    message(FATAL_ERROR "MULTIBLEND_BUILD_EXE requires jpeg, png and tiff support.")
  endif()
endif()


include(CTest)

set(CMAKE_CXX_STANDARD 20)

add_library(MultiblendLib
  "src/aligned_ptr.cpp"
  "src/functions.cpp"
  "src/image.cpp"
  "src/mapalloc.cpp"
  "src/multiblend.cpp"
  "src/pnger.cpp"
  "src/pyramid.cpp"
  "src/threadpool.cpp"
)

target_include_directories(MultiblendLib PUBLIC 
  "include"
)

find_package(spdlog REQUIRED)
target_link_libraries(MultiblendLib spdlog::spdlog)

if(MULTIBLEND_LIB_WITH_JPEG)
  find_package(JPEG REQUIRED)
  target_link_libraries(MultiblendLib JPEG::JPEG)
  target_compile_definitions(MultiblendLib PUBLIC MULTIBLEND_WITH_JPEG)
endif()

if(MULTIBLEND_LIB_WITH_PNG)
  find_package(PNG REQUIRED)
  target_link_libraries(MultiblendLib PNG::PNG)
  target_compile_definitions(MultiblendLib PUBLIC MULTIBLEND_WITH_PNG)
endif()

if(MULTIBLEND_LIB_WITH_TIFF)
  find_package(TIFF REQUIRED)
  target_link_libraries(MultiblendLib TIFF::TIFF)
  target_compile_definitions(MultiblendLib PUBLIC MULTIBLEND_WITH_TIFF)
  target_sources(MultiblendLib PRIVATE "src/tiff.cpp")
endif()

if(MULTIBLEND_BUILD_EXE)
  add_executable(Multiblend
    "src/main.cpp"
  )

  find_package(JPEG REQUIRED)
  find_package(PNG REQUIRED)
  find_package(TIFF REQUIRED)

  target_link_libraries(Multiblend
    JPEG::JPEG
    PNG::PNG
    TIFF::TIFF
    MultiblendLib
  )
endif()

if(NOT MSVC)
    target_compile_options(MultiblendLib PRIVATE -msse4.1)
endif()

if(MULTIBLEND_WITH_ASAN)
  set(ASAN_CXX_OPTIONS "-fsanitize=address;-g")
  set(ASAN_LINK_OPTIONS "-fsanitize=address")

  target_compile_options(MultiblendLib PUBLIC "${ASAN_CXX_OPTIONS}")
  target_link_options(MultiblendLib PUBLIC "${ASAN_LINK_OPTIONS}")

  if(MULTIBLEND_BUILD_EXE)
    target_compile_options(Multiblend PUBLIC "${ASAN_CXX_OPTIONS}")
    target_link_options(Multiblend PUBLIC "${ASAN_LINK_OPTIONS}")
  endif()
endif()

if(BUILD_TESTING AND ${MULTIBLEND_ROOT_PROJECT})
  add_subdirectory(tests)
endif()
