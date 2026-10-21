cmake_minimum_required(VERSION 3.18)
project(timberwolf C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# ---- Interface library ----
add_library(timberwolf INTERFACE)
target_include_directories(timberwolf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# ---- CHEAPDAN (auto-fetches FFTW3 + CHEAP dependency) ----
include(FetchContent)

FetchContent_Declare(
    cheapdan
    GIT_REPOSITORY https://github.com/kevinfling/CHEAPDAN.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(cheapdan)

target_link_libraries(timberwolf INTERFACE cheapdan)

# ---- Subdirectories ----
enable_testing()
add_subdirectory(tests)
add_subdirectory(benchmarks)
add_subdirectory(examples)

# ---- Install ----
# CHEAPDAN is a fetched dependency and not included in the export set;
# link it privately at install time via an explicit find_dependency call.
install(TARGETS timberwolf EXPORT timberwolfTargets)
install(FILES include/timberwolf.h DESTINATION include)
