cmake_minimum_required(VERSION 3.21)

find_package (Python COMPONENTS Interpreter REQUIRED)

add_test(NAME MultiblendTest COMMAND 
  ${Python_EXECUTABLE} -m unittest discover -s ${CMAKE_CURRENT_SOURCE_DIR} -v)

execute_process(COMMAND 
  ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data)
