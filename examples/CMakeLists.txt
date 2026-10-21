cmake_minimum_required(VERSION 3.18)

add_executable(basic_health_index basic_health_index.c)
target_link_libraries(basic_health_index PRIVATE timberwolf)
