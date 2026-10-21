cmake_minimum_required(VERSION 3.18)

add_executable(bench_process bench_process.c)
target_link_libraries(bench_process PRIVATE timberwolf)
