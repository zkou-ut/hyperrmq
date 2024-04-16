add_executable(benchmark_arithmetic_decoding benchmark_arithmetic_decoding.cpp)
target_link_libraries(benchmark_arithmetic_decoding hyperrmq)

add_executable(change_n change_n.cpp)
target_link_libraries(change_n hyperrmq ferradarmq sdsl)

add_executable(change_width change_width.cpp)
target_link_libraries(change_width hyperrmq ferradarmq sdsl)

add_executable(lcp lcp.cpp)
target_link_libraries(lcp hyperrmq ferradarmq sdsl)
