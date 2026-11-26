cmake_minimum_required(VERSION 3.15)
project(matmpc_core LANGUAGES C)

add_library(matmpc_core STATIC erk_kernel.c irk_kernel.c irk_dae_kernel.c recover_kernel.c matmpc_la.c matmpc_lu.c
    condense_affine.c condense_matrix.c condense_mb.c qp_gauss_newton.c)
target_include_directories(matmpc_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(matmpc_core PUBLIC c_std_99)
set_target_properties(matmpc_core PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
if(NOT MSVC)
    target_link_libraries(matmpc_core PUBLIC m)
endif()

include(CTest)
if(BUILD_TESTING)
    foreach(test_name IN ITEMS recover erk_kernel irk_kernel irk_dae_kernel matmpc_la matmpc_lu
        condense_affine condense_matrix condense_mb qp_gauss_newton qp_implicit)
        add_executable(test_${test_name} "../tests/native/test_${test_name}.c")
        target_link_libraries(test_${test_name} PRIVATE matmpc_core)
        set_target_properties(test_${test_name} PROPERTIES C_EXTENSIONS OFF)
        add_test(NAME ${test_name} COMMAND test_${test_name})
    endforeach()
endif()
