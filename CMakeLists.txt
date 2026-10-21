cmake_minimum_required(VERSION 3.18)
project(timberwolf LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(TIMBERWOLF_BUILD_TESTS "Build tests" ON)
option(TIMBERWOLF_BUILD_EXAMPLES "Build examples" ON)
option(TIMBERWOLF_BUILD_BENCHMARKS "Build benchmarks" ON)

add_library(timberwolf INTERFACE)
target_include_directories(timberwolf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(FetchContent)
# Legacy dependency: pin a specific commit for reproducible releases.
FetchContent_Declare(cheapdan
    GIT_REPOSITORY https://github.com/kevinfling/CHEAPDAN.git
    GIT_TAG master
)
FetchContent_MakeAvailable(cheapdan)
target_link_libraries(timberwolf INTERFACE cheapdan)

if(TIMBERWOLF_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
if(TIMBERWOLF_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
if(TIMBERWOLF_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

install(TARGETS timberwolf EXPORT timberwolfTargets)
install(FILES include/timberwolf.h DESTINATION include)
