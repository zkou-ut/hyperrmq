cmake_minimum_required(VERSION 3.13)
project(hyperrmq CXX)

set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wno-sign-compare -Wextra -Wno-type-limits")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")

# hyperrmq
add_subdirectory(src)

# Ferrada & Navarro RMQ library
add_library(ferradarmq SHARED rmq/includes/Basic_rmq.cpp rmq/RMQRMM64.cpp)
target_include_directories(ferradarmq PUBLIC rmq/includes)

# sdsl
add_subdirectory(sdsl-lite)

# manually link include directories for sub libraries that
target_include_directories(gtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/sdsl-lite/external/googletest/include)
target_include_directories(gtest_main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/sdsl-lite/external/googletest/include)
target_include_directories(divsufsort PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/external/libdivsufsort/include)
target_include_directories(divsufsort64 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/external/libdivsufsort/include)


add_subdirectory(experiment)
add_subdirectory(test)

add_executable(simple-test simple_test.cpp)
target_link_libraries(simple-test hyperrmq ferradarmq sdsl)
