cmake_minimum_required(VERSION 2.8)
project(whiteboard)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11" CACHE STRING "" FORCE)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

set(SRC_PATH ${CMAKE_SOURCE_DIR}/src)
include_directories(${SRC_PATH})
include_directories(/usr/local/include/)

file(GLOB lib_src ${SRC_PATH}/*.cpp)
add_library(whiteboard ${lib_src})
add_executable(main systest/main.cpp)
target_link_libraries(main whiteboard ${OpenCV_LIBS})
