file(GLOB SOURCES CONFIGURE_DEPENDS *.cpp)

add_library(hyperrmq SHARED ${SOURCES})

target_include_directories(hyperrmq PUBLIC ${PROJECT_SOURCE_DIR}/include)
