set(PROGRAM hyperrmq-test)

file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/*.cpp)

add_executable(${PROGRAM} ${TEST_SOURCES})

target_link_libraries(${PROGRAM} PRIVATE hyperrmq gtest)
