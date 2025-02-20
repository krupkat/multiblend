cmake_minimum_required(VERSION 3.21)

execute_process(COMMAND 
  ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data)

if(MULTIBLEND_BUILD_EXE)
  find_package (Python3 COMPONENTS Interpreter REQUIRED)

  add_test(NAME MultiblendTest COMMAND 
    ${Python3_EXECUTABLE} -m unittest discover -s ${CMAKE_CURRENT_SOURCE_DIR} -v)
endif()

find_package(Catch2 REQUIRED)
find_package(cereal REQUIRED)

add_executable(MultiblendLibTest
  multiblend_test.cpp
)

target_link_libraries(MultiblendLibTest
  PRIVATE
    MultiblendLib
    Catch2::Catch2WithMain
    cereal::cereal
)

include(Catch)
catch_discover_tests(MultiblendLibTest
  WORKING_DIRECTORY "$<TARGET_FILE_DIR:MultiblendLibTest>"
)
