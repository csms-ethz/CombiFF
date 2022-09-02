cmake_minimum_required (VERSION 2.8.12)

set(shared_path ${CMAKE_CURRENT_SOURCE_DIR}/shared)

#to choose Debug build type, in build, use cmake -DCMAKE_BUILD_TYPE="Debug" ..
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release. If you want to use debug mode, use -DCMAKE_BUILD_TYPE=\"Debug\"")
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message("Setting CMAKE_INSTALL_PREFIX to ${CMAKE_SOURCE_DIR} instead of default ${CMAKE_INSTALL_PREFIX}. You can change the path by using -DCMAKE_INSTALL_PREFIX=<path>")
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "install in CMake source directory" FORCE)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wpedantic -O3 -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wpedantic -Wconversion -O3 -std=c++11 -g")

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(shared)
add_subdirectory(enu)
add_subdirectory(cnv)
add_subdirectory(tbl)
