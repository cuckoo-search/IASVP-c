cmake_minimum_required(VERSION 3.3)

project(cuckoo_search_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -march=native -mtune=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -flto -ffast-math -funroll-all-loops -march=native -mtune=native")

set(SOURCE_FILES main.cpp)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CMAKE_CXX_COMPILER g++-5)

include_directories(includes IASVP)

add_executable(cuckoo_search_cpp ${SOURCE_FILES})

target_link_libraries(cuckoo_search_cpp m lapack cblas blas)
