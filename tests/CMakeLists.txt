cmake_minimum_required(VERSION 3.18)

include(FetchContent)

FetchContent_Declare(
    unity
    GIT_REPOSITORY https://github.com/ThrowTheSwitch/Unity.git
    GIT_TAG        v2.6.0
)
FetchContent_MakeAvailable(unity)

# Unity must be compiled with UNITY_INCLUDE_DOUBLE so UnityAssertDoublesWithin is linked in
target_compile_definitions(unity PUBLIC UNITY_INCLUDE_DOUBLE)

function(add_tw_test name src)
    add_executable(${name} ${src})
    target_compile_definitions(${name} PRIVATE UNITY_INCLUDE_DOUBLE)
    target_link_libraries(${name} PRIVATE timberwolf unity)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_tw_test(test_robust_stats  test_robust_stats.c)
add_tw_test(test_wolf_ewma     test_wolf_ewma.c)
add_tw_test(test_pipeline      test_pipeline.c)
