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


add_subdirectory(experiment)
