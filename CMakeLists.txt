cmake_minimum_required(VERSION 3.16)
project(TARAL_LP LANGUAGES CXX)

# Source order below matches the previous CMakeLists exactly (link order affects the binary).
# Optional switches. Both default to OFF, so a plain configure builds exactly the same
# `taral` target with the same flags as before.
option(TARAL_WARNINGS   "Build with -Wall -Wextra -Wpedantic" OFF)
option(TARAL_BUILD_TESTS "Build the C++ unit tests and register them with ctest" OFF)

set(TARAL_ENGINE_SOURCES
    src/crossover.cpp
    src/dual.cpp
    src/explain.cpp
    src/ipm.cpp
    src/kkt_gate.cpp
    src/lu.cpp
    src/main.cpp
    src/milp.cpp
    src/mps.cpp
    src/nonoptimal_certificate.cpp
    src/presolve.cpp
    src/simplex.cpp
)

add_executable(taral ${TARAL_ENGINE_SOURCES})
target_compile_features(taral PRIVATE cxx_std_17)
set_target_properties(taral PROPERTIES CXX_EXTENSIONS OFF)
if(TARAL_WARNINGS)
    target_compile_options(taral PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(TARAL_BUILD_TESTS)
    enable_testing()

    set(TARAL_LIB_SOURCES ${TARAL_ENGINE_SOURCES})
    list(REMOVE_ITEM TARAL_LIB_SOURCES src/main.cpp)

    # Engine sources (everything except main.cpp) compiled once and shared by the unit tests.
    # The tests use assert(), so NDEBUG is undefined for them even in Release builds.
    add_library(taral_core OBJECT ${TARAL_LIB_SOURCES})
    target_compile_features(taral_core PUBLIC cxx_std_17)
    set_target_properties(taral_core PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(taral_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
    target_compile_options(taral_core PUBLIC -UNDEBUG)

    function(taral_unit_test name source)
        add_executable(${name} ${source} $<TARGET_OBJECTS:taral_core>)
        target_link_libraries(${name} PRIVATE taral_core)
        add_test(NAME ${name} COMMAND ${name} ${ARGN})
    endfunction()

    taral_unit_test(kkt_gate_tests tests/kkt_gate_tests.cpp)
    taral_unit_test(certificate_quality_api_tests benchmarks/certificate_quality_api_tests.cpp)
    taral_unit_test(infeasibility_explanation_api_tests benchmarks/infeasibility_explanation_api_tests.cpp)
    taral_unit_test(nonoptimal_certificate_api_tests benchmarks/nonoptimal_certificate_api_tests.cpp)
    taral_unit_test(presolve_api_tests benchmarks/presolve_api_tests.cpp)
    taral_unit_test(deterministic_solve_api_tests benchmarks/deterministic_solve_api_tests.cpp
                    ${CMAKE_SOURCE_DIR}/tests/fixtures/presolve_fixedbig.mps)

    # CLI and regression scripts; they run from the repository root.
    foreach(script cli_test parse_deadline_test presolve_work_limit_test)
        add_test(NAME ${script}
                 COMMAND sh ${CMAKE_SOURCE_DIR}/tests/${script}.sh $<TARGET_FILE:taral>
                 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endforeach()
    add_test(NAME ipm_gate
             COMMAND bash ${CMAKE_SOURCE_DIR}/tests/ipm_gate/check.sh $<TARGET_FILE:taral>
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
