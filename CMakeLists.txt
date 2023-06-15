cmake_minimum_required(VERSION 3.21)
project(Multiblend)

cmake_host_system_information(RESULT MULTIBLEND_HAS_SSE QUERY HAS_SSE)

set(MULTIBLEND_ROOT_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MULTIBLEND_ROOT_PROJECT ON)
endif()

include (CMakeDependentOption)

OPTION(MULTIBLEND_WITH_ASAN "Build with address sanitizer" OFF)
OPTION(MULTIBLEND_WITH_TSAN "Build with thread sanitizer" OFF)
OPTION(MULTIBLEND_BUILD_EXE "Build executable" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_LIB_WITH_JPEG "Build with JPEG support" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_LIB_WITH_PNG "Build with PNG support" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_LIB_WITH_TIFF "Build with TIFF support" ${MULTIBLEND_ROOT_PROJECT})
OPTION(MULTIBLEND_ARM_OPTIMIZED "Enable optimizations for ARM" OFF)

set(MULTIBLEND_ARM_ARCHITECTURE "armv8-a+simd"
  CACHE STRING "ARM architecture target"
)

cmake_dependent_option(MULTIBLEND_WITH_OPENMP_SIMD "Build SIMDe with OpenMP" 
  ON "MULTIBLEND_ARM_OPTIMIZED" OFF
)

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
  "src/logging.cpp"
  "src/mapalloc.cpp"
  "src/multiblend.cpp"
  "src/pnger.cpp"
  "src/pyramid.cpp"
)

target_include_directories(MultiblendLib PUBLIC 
  "include"
)

find_package(spdlog REQUIRED)
target_link_libraries(MultiblendLib PUBLIC spdlog::spdlog)

find_path(SIMDE_INCLUDE_DIRS "simde/arm/neon.h")
target_include_directories(MultiblendLib PRIVATE ${SIMDE_INCLUDE_DIRS})

if(MULTIBLEND_LIB_WITH_JPEG)
  find_package(JPEG REQUIRED)
  target_link_libraries(MultiblendLib PUBLIC JPEG::JPEG)
  target_compile_definitions(MultiblendLib PUBLIC MULTIBLEND_WITH_JPEG)
endif()

if(MULTIBLEND_LIB_WITH_PNG)
  find_package(PNG REQUIRED)
  target_link_libraries(MultiblendLib PUBLIC PNG::PNG)
  target_compile_definitions(MultiblendLib PUBLIC MULTIBLEND_WITH_PNG)
endif()

if(MULTIBLEND_LIB_WITH_TIFF)
  find_package(TIFF REQUIRED)
  target_link_libraries(MultiblendLib PUBLIC TIFF::TIFF)
  target_compile_definitions(MultiblendLib PUBLIC MULTIBLEND_WITH_TIFF)
  target_sources(MultiblendLib PRIVATE "src/tiff.cpp")
endif()

if(MULTIBLEND_BUILD_EXE)
  add_executable(Multiblend
    "src/main.cpp"
  )

  target_link_libraries(Multiblend PRIVATE
    MultiblendLib
  )
endif()

if(MULTIBLEND_ARM_OPTIMIZED)
  target_compile_options(MultiblendLib PRIVATE "-march=${MULTIBLEND_ARM_ARCHITECTURE}")
elseif(MULTIBLEND_HAS_SSE AND NOT MSVC)
  target_compile_options(MultiblendLib PRIVATE -msse4.1)
endif()

if(MULTIBLEND_WITH_OPENMP_SIMD)
  target_compile_options(MultiblendLib PRIVATE "-fopenmp-simd")
  target_compile_definitions(MultiblendLib PRIVATE SIMDE_ENABLE_OPENMP)
endif()

function(multiblend_build_with_sanitizer sanitizer)
  set(SANITIZER_CXX_OPTIONS "-fsanitize=${sanitizer};-g")
  set(SANITIZER_LINK_OPTIONS "-fsanitize=${sanitizer}")

  target_compile_options(MultiblendLib PUBLIC "${SANITIZER_CXX_OPTIONS}")
  target_link_options(MultiblendLib PUBLIC "${SANITIZER_LINK_OPTIONS}")

  if(MULTIBLEND_BUILD_EXE)
    target_compile_options(Multiblend PUBLIC "${SANITIZER_CXX_OPTIONS}")
    target_link_options(Multiblend PUBLIC "${SANITIZER_LINK_OPTIONS}")
  endif()
endfunction()

if(MULTIBLEND_WITH_ASAN)
  multiblend_build_with_sanitizer(address)
endif()

if(MULTIBLEND_WITH_TSAN)
  multiblend_build_with_sanitizer(thread)
endif()

if(BUILD_TESTING AND ${MULTIBLEND_ROOT_PROJECT})
  add_subdirectory(tests)
endif()
